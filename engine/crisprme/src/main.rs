use clap::{Parser, Subcommand};
use crisprme_core::common::alignment::{Alignment, visualize};
use crisprme_core::common::guide::Guide;
use crisprme_core::common::iupac::Iupac;
use crisprme_core::common::sequence::Sequence;
use crisprme_core::common::cigarx::{Cigarx, CigarxOp};
use crisprme_core::memory::arena::Arena;
use crisprme_core::pipeline::engine::hybrid::HybridEngine;
use crisprme_core::pipeline::PipelineDescriptor;
use crisprme_core::storage::reader::BinaryPositionReader;
use crisprme_core::storage::reader::BinarySequenceReader;
use rayon::prelude::*;
use serde_json::Value;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use std::path::PathBuf;
use tracing::error;
use tracing::info;

use crisprme_core::utils;

// Use custom small allocator
use mimalloc::MiMalloc;
use ahash::AHashMap;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

type JsonMap = serde_json::Map<String, Value>;

mod common;
use common::*;

#[derive(Parser)]
#[command(version)]
struct Cli {
    /// Show logs
    #[arg(short, long)]
    verbose: bool,

    /*
    /// Threads configuration
    #[command(flatten)]
    threads_config: ThreadsConfig,
    */

    /// Memory configuration
    #[command(flatten)]
    memory_config: MemoryConfig,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Show alignment of two sequences with a CIGARX
    Cigarx {
        guide: String,
        target: String,
        cigarx: String,
        offset: usize,
    },

    /// Preprocess a FASTA file to obtain binary packed sequences and positions
    Preprocess {
        /// Filepath to the input fasta file
        input: PathBuf,
        /// Length of the sequences
        sequence_len: usize,
        /// Delta between sequences
        delta: usize,
    },

    /// Preprocess a ids + sequence file to obtain binary packed sequences and positions
    PreprocessList {
        /// Filepath to the input fasta file
        input: PathBuf,
        /// Length of the sequences
        sequence_len: usize,
    },

    /// Show all stored sequences and ids inside a split binary file
    Sequences {
        /// Name of the target dataset
        input: String,
        /// Length of the sequences
        sequence_len: usize,
    },

    Mine {
        /// Name of the target dataset
        input: String,
        /// Length of the sequences
        sequence_len: usize,
        /// Guide sequence to align
        guide: String,
        /// Thresholds configuration
        #[command(flatten)]
        thresholds: CliThresholds,
        /// Filepath to output alignment file
        output: PathBuf,
    },

    Alignments {
        /// Filepath to mined binary file
        input: PathBuf,

        /// Show only positions
        #[clap(long, short, action)]
        positions: bool
    },

    Results {
        /// Name of the target dataset
        input: String,
        /// Filepath of the alignments binary file
        alignments: PathBuf,
        /// Length of the sequences
        sequence_len: usize,
        /// Guide used for the alignment
        guide: String,
        /// Skip wildcards
        #[clap(long, short, action)]
        skip_wildcards: bool,
        /// Show the some complete alignments
        #[clap(long, action)]
        preview: bool
    }
}

fn main() {
    let cli = Cli::parse();

    tracing_subscriber::fmt()
        .compact()
        .with_target(false)
        .with_thread_ids(true)
        .with_max_level(tracing::Level::TRACE)
        .init();

    match &cli.command {
        Commands::Cigarx {
            guide,
            target,
            cigarx,
            offset,
        } => {
            utils::visualize::cigar(
                guide.as_bytes(),
                target.as_bytes(),
                cigarx.as_bytes(),
                *offset,
            );
        }

        Commands::Preprocess {
            input,
            sequence_len,
            delta,
        } => {
            // Create output folder
            let stem = input.file_stem().unwrap();
            let output_folder = format!("./data/ready/{}", stem.to_string_lossy());
            std::fs::create_dir_all(&output_folder)
                .expect("unable to create output folder structure");

            // Create witers
            let mut pos_writer = BufWriter::new(
                File::create(format!("{}/positions.bin", &output_folder))
                    .expect("unable to create binary position file"),
            );
            let mut seq_writer = BufWriter::new(
                File::create(format!("{}/sequences.bin", &output_folder))
                    .expect("unable to create binary sequence file"),
            );

            // Load entire fasta file
            let fasta = std::fs::read(input).expect("unable to read FASTA file");

            // Skip header ">chr..."
            let initial = fasta.iter().position(|&e| e == b'\n').unwrap_or_default();

            // Convert to IUPAC and remove newlines
            let fasta: Vec<Iupac> = fasta
                .into_iter()
                .skip(initial)
                .filter(|&e| e != b'\n')
                .map(Iupac::from_ascii)
                .collect();

            // Find first occurence of a non N character
            let beg = match fasta.iter().skip(10).position(|e| !e.is_wildcard()) {
                Some(start) => { 
                    info!("First non N character at {}", start);
                    start 
                },
                None => {
                    error!("FASTA file has only Ns");
                    return;
                }
            };

            // Find the length of the valid characters
            // SAFETY: If we have a start index then this cannot fail
            let end = fasta.iter().rposition(|e| !e.is_wildcard()).unwrap() + 1;
            info!("Ending N character at {}", end);

            // Work only on valid region
            let fasta = &fasta[beg..end];
            let mut iter = fasta.windows(*sequence_len).step_by(*delta).enumerate();
            let chunks: Vec<Vec<_>> = {
                let mut output = Vec::new();
                loop {
                    let chunk: Vec<_> = iter.by_ref().take(1_000_000).collect();
                    if chunk.is_empty() {
                        break;
                    }
                    output.push(chunk);
                }
                output
            };

            info!("generated {} chunks for parallel processing", chunks.len());

            // Process chunks in parallel
            let sequences = chunks
                .iter()
                .par_bridge()
                .map(|chunk| {
                    let mut sequences: BTreeMap<&[Iupac], Vec<u32>> = BTreeMap::new();
                    for (i, seq) in chunk {
                        // Add position to sequence key
                        let pos = (beg + i * (*delta)) as u32;
                        sequences.entry(seq).or_default().push(pos);
                    }
                    sequences
                })
                .reduce(BTreeMap::new, |mut a, b| {
                    for (k, v) in b {
                        a.entry(k).and_modify(|s| s.extend(v.clone())).or_insert(v);
                    }
                    a
                });

            let mut sequence_count = sequences.keys().len();
            info!("found {} unique sequences", sequence_count);

            // Write to files
            for (seq, positions) in sequences.iter() {
                sequence_count -= 1;

                // Write ids: [len:u32][pos[]:[u32]]
                let len = positions.len() as u32;
                pos_writer.write_all(&len.to_le_bytes()).unwrap();
                for pos in positions {
                    pos_writer.write_all(&pos.to_le_bytes()).unwrap();
                }

                // Write sequences
                let bytes =
                    unsafe { std::slice::from_raw_parts(seq.as_ptr() as *const u8, seq.len()) };
                seq_writer.write_all(bytes).unwrap();
            }

            assert!(sequence_count == 0);
        },

        Commands::PreprocessList { input, sequence_len } => {

            // Open input file
            let reader = BufReader::new(File::open(&input)
                .expect("unable to open file"));
            
            // Create output folder
            let stem = input.file_stem().unwrap();
            let output_folder = format!("./data/ready/{}", stem.to_string_lossy());
            std::fs::create_dir_all(&output_folder)
                .expect("unable to create output folder structure");

            // Create witers
            let mut pos_writer = BufWriter::new(
                File::create(format!("{}/positions.bin", &output_folder))
                    .expect("unable to create binary position file"),
            );
            let mut seq_writer = BufWriter::new(
                File::create(format!("{}/sequences.bin", &output_folder))
                    .expect("unable to create binary sequence file"),
            );

            info!("reading dataset {} into memory", stem.to_string_lossy());

            // NOTE: We assume an ASCII encoded file
            let dataset: Vec<u8> = std::fs::read(input).expect("unable to read list file");
            let lines: Vec<&[u8]> = dataset.split(|&b| b == b'\n').collect();

            // Create the maximum chunk size for the number of threads
            let parallelism = rayon::current_num_threads();
            let chunk_size = (lines.len() + parallelism - 1) / parallelism;

            // NOTE: This uses a faster hashing function
            let mut result: AHashMap<&[u8], Vec<u32>> = AHashMap::new();

            info!("aggregating across {} threads", parallelism);

            let chunk_maps: Vec<(BTreeMap<_, _>, u32)> = lines.par_chunks(chunk_size)
                .map(|chunk| {

                    let mut result: BTreeMap<&[u8], Vec<u32>> = BTreeMap::new();
                    let mut lines_count = 0;
                    for line in chunk {
                        if let Some(mid) = line.iter().position(|&b| b == b'\t') {
                            lines_count += 1;

                            let position = &line[..mid];
                            let position = unsafe {
                                std::str::from_utf8_unchecked(position)
                                    .parse::<u32>().unwrap()
                            };

                            let sequence_bytes = &line[mid+1..];
                            result.entry(sequence_bytes).or_default()
                                .push(position);
                        }
                    }

                    info!("completed partial aggregation");
                    (result, lines_count)
                })
                .collect();
            
            info!("aggregating partials");

            // Aggregate the map from each thread
            // NOTE: Global allocation is done in a single thread to reduce memory usage
            let (result, total_lines) = chunk_maps.into_iter()
                .reduce(
                    |(mut acc, lines), (chunk_map, chunk_lines)| {
                        for (k, mut v) in chunk_map {
                            acc.entry(k)
                                .and_modify(|r| r.append(&mut v))
                                .or_insert(v);
                        }
                        (acc, lines + chunk_lines)
                    }
                ).unwrap();

            info!("found {} unique sequences compared to the initial {}", result.len(), total_lines);
            info!("writing packed data");

            let mut positions_count = 0;
            let mut sequence_min_positions = usize::MAX;
            let mut sequence_max_positions = usize::MIN;

            // Write sequences and ids to file
            let mut scratch: [Iupac; 128] = [Iupac::from_ascii(b'N'); 128];
            for (seq, positions) in result.iter() {
                
                sequence_min_positions = sequence_min_positions.min(positions.len());
                sequence_max_positions = sequence_max_positions.max(positions.len());

                // Write ids: [len:u32][pos[]:[u32]]
                let len = positions.len() as u32;
                pos_writer.write_all(&len.to_le_bytes()).unwrap();
                for pos in positions {
                    pos_writer.write_all(&pos.to_le_bytes()).unwrap();
                    positions_count += 1;
                }

                // Convert the seq from u8 to Iupac
                for (i, c) in seq.iter().enumerate() {
                    scratch[i] = Iupac::from_ascii(*c);
                }

                // Write sequences
                // SAFETY: Iupac is repr(C) and is a u8
                let bytes =
                    unsafe { std::slice::from_raw_parts(scratch.as_ptr() as *const u8, *sequence_len) };

                seq_writer.write_all(bytes).unwrap();
            }

            info!("{} written records", positions_count);
            info!("max positions for a sequence: {}", sequence_max_positions);
            info!("min positions for a sequence: {}", sequence_min_positions);
            assert_eq!(positions_count, total_lines);

            /*
            let slice = &fasta[beg..end];
            let mut iter = slice.windows(*sequence_len).step_by(*delta).enumerate();
            let chunks: Vec<Vec<_>> = {
                let mut output = Vec::new();
                loop {
                    let chunk: Vec<_> = iter.by_ref().take(1_000_000).collect();
                    if chunk.is_empty() {
                        break;
                    }
                    output.push(chunk);
                }
                output
            };


            let mut sequences: BTreeMap<Vec<Iupac>, Vec<u32>> = BTreeMap::new();
            for line in reader.lines() {
                let line = line.unwrap();

                let parts: Vec<&str> = line.split('\t').collect();
                let id = parts[0].parse::<u32>().unwrap();
                let sequence: Vec<Iupac> = parts[1].chars()
                    .map(|e| e as u8)
                    .map(Iupac::from_ascii)
                    .collect();

                assert_eq!(sequence.len(), *sequence_len);
                //println!("id: {}, sequence: {}", id, Sequence::new(&sequence));
                sequences.entry(sequence)
                    .or_insert_with(Vec::new)
                    .push(id);
            }

            // Write sequences and ids to file
            let mut sequence_count = 0;
            for (seq, positions) in sequences.iter() {
                sequence_count -= 1;

                // Write ids: [len:u32][pos[]:[u32]]
                let len = positions.len() as u32;
                pos_writer.write_all(&len.to_le_bytes()).unwrap();
                for pos in positions {
                    pos_writer.write_all(&pos.to_le_bytes()).unwrap();
                }

                // Write sequences
                let bytes =
                    unsafe { std::slice::from_raw_parts(seq.as_ptr() as *const u8, seq.len()) };
                seq_writer.write_all(bytes).unwrap();
            }

            info!("stored {} sequences", sequence_count);
            */
        },

        Commands::Mine { input, output, guide, sequence_len, thresholds } => {

            let input_path = format!("./data/ready/{}/sequences.bin", input);
            let input = Path::new(&input_path);

            let engine = HybridEngine {};
            engine.execute(PipelineDescriptor {
                guide: Guide::from(guide.as_str()),
                sequence_len: *sequence_len,
                sequence_batch_size: cli.memory_config.sequence_batch_size,
                alignment_batch_size: cli.memory_config.alignment_batch_size,
                sequence_file: input.to_path_buf(),
                output_file: output.clone(),
                thresholds: thresholds.into(),
                mutation_max: 100,
            });
        }

        Commands::Sequences { input, sequence_len } => {
            let pos_file = File::open(format!("./data/ready/{}/positions.bin", input))
                .expect("unable to open ids file");

            let seq_file = File::open(format!("./data/ready/{}/sequences.bin", input))
                .expect("unable to open sequence file");

            let mut arena = Arena::alloc(1024 * 1024 * 1024);
            arena.scoped(|mem| {
                let pos_reader = BinaryPositionReader::new(&mem, BufReader::new(pos_file));
                let seq_reader = BinarySequenceReader::new(&mem, BufReader::new(seq_file), *sequence_len);

                for (i, (ids, seq)) in pos_reader.zip(seq_reader).enumerate() {
                    let seq = Sequence::new(seq);
                    println!("index {}:", i);
                    println!("   sequence (IUPAC): {}", seq);
                    println!("   positions  (u32): {:?}", ids);
                }
            });
        }

        Commands::Alignments { input, positions } => {
            let mut reader = BufReader::new(File::open(input).expect("file not found"));

            let mut positions_set = HashSet::new();
            let mut buffer = vec![0u8; size_of::<Alignment>()];
            loop {
                match reader.read_exact(&mut buffer) {
                    Ok(()) => {}
                    Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                    Err(e) => println!("ERR: {e}"),
                };

                let alignment = unsafe { &*(buffer.as_ptr() as *const Alignment) };
                if !*positions {
                    println!("{alignment:?}");
                } else {
                    let global_pos = alignment.id + alignment.offset as u32;
                    positions_set.insert(global_pos);
                }
            }

            if *positions {
                for p in positions_set {
                    println!("{p}");
                }
            }
        },

        Commands::Results { input, alignments, sequence_len, guide, preview, skip_wildcards } => {
           
            let positive_guide = Guide::from(guide.as_str());
            let negative_guide = positive_guide.reverse_complement();

            let mut alig_reader = BufReader::new(File::open(alignments)
                .expect("file not found"));

            let seq_file = File::open(format!("./data/ready/{}/sequences.bin", input))
                .expect("unable to open sequence file");

            let pos_file = File::open(format!("./data/ready/{}/positions.bin", input))
                .expect("unable to open ids file");

            let mut sequence_ids = HashMap::new();
            let mut buffer = vec![0u8; size_of::<Alignment>()];
            loop {
                match alig_reader.read_exact(&mut buffer) {
                    Ok(()) => {}
                    Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                    Err(e) => println!("ERR: {e}"),
                };

                // SAFETY: We must clone! Otherwise we keep only a reference to the last alignment
                // for all values of the hashmap!
                let alignment = unsafe { &*(buffer.as_ptr() as *const Alignment) };
                sequence_ids.entry(alignment.id as u32)
                    .or_insert_with(Vec::new)
                    .push(alignment.clone());
            }

            /*
            for (k, v) in sequence_ids {
                println!("id: {k}, align: {v:?}");
            }
            */

            //let positive_guide_bytes = guide.bytes();
            //let negative_guide_bytes = None;

            let mut arena = Arena::alloc(1024 * 1024 * 1024);
            arena.scoped(|mem| {
                let pos_reader = BinaryPositionReader::new(&mem, BufReader::new(pos_file));
                let seq_reader = BinarySequenceReader::new(&mem, BufReader::new(seq_file), *sequence_len);

                println!("position,offset,strand,cigarx,guide,sequence");
                for (i, (positions, seq)) in pos_reader.zip(seq_reader).enumerate() {
                    let seq_id = i as u32;
                    if let Some(alignments) = sequence_ids.get(&seq_id) {
                        let sequence = Sequence::new(seq);

                        // Skip all N results if requested
                        if *skip_wildcards {
                            if sequence.mutation_score() == (4 * sequence_len) as u32 {
                                continue;
                            }
                        }

                        // Print all alignment at all positions
                        // POS, OFFSET, STRAND, CIGARX, SEQ
                        for align in alignments {
                            for pos in &positions {

                                let guide = if align.strand == b'+' {
                                    &positive_guide
                                } else {
                                    &negative_guide
                                };

                                // Print complete solution
                                let mut qline = String::new();
                                let mut cline = String::new();
                                let mut tline = String::new();
                                
                                let mut qidx: usize = 0;
                                let mut tidx: usize = 0;

                                // Add prefix
                                for i in 0..align.offset {
                                    tline.push(sequence[tidx].to_utf8());
                                    tidx += 1;
                                }

                                // Alignments
                                for op in align.cigarx.operations() {
                                    cline.push(op.to_utf8());
                                    match op {
                                        CigarxOp::Match | CigarxOp::Mismatch => {
                                            tline.push(sequence[tidx].to_utf8());
                                            qline.push(guide[qidx].to_utf8());
                                            qidx += 1;
                                            tidx += 1;
                                        },
                                        CigarxOp::Deletion => {
                                            qline.push(guide[qidx].to_utf8());
                                            tline.push('-');
                                            qidx += 1;
                                        },
                                        CigarxOp::Insertion => {
                                            tline.push(sequence[tidx].to_utf8());
                                            qline.push('-');
                                            tidx += 1;
                                        }
                                    }
                                }
                                
                                // Add target suffix (unaligned)
                                while tidx < *sequence_len {
                                    tline.push(sequence[tidx].to_utf8());
                                    tidx += 1;
                                }

                                println!("{},{},{},{},{},{}",
                                    pos,
                                    align.offset,
                                    align.strand as char,
                                    align.cigarx,
                                    qline,
                                    tline
                                );

                                /*
                                println!("{},{},{},{},{}", 
                                    pos,
                                    align.offset,
                                    align.strand as char,
                                    align.cigarx,
                                    seq
                                );
                                */
                            }
                        }
                    }
                }
            });
        }
    }
}
