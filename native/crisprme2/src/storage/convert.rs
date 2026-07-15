use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::sequence::iupac::Iupac;

/// Read a text genome file and convert it to a packed binary file
pub fn binary_block_pack_soa<P: AsRef<Path>>(input: P, output: P, seq_len: usize) {
    let input = input.as_ref();
    assert_eq!("txt", input.extension().unwrap());

    let filename = input.file_name().unwrap();
    let filename = filename.to_str().unwrap();
    let (basename, _) = filename.split_at(filename.find('.').unwrap());

    let output = output.as_ref().join(basename);

    let output_file = std::fs::File::create(&output).expect("unable to create output file");
    let mut writer = BufWriter::new(output_file);

    let mut total_records = 0;

    // Read all ids
    let reader = BufReader::new(File::open(input).expect("unable to open file"));
    for content in reader.lines() {
        let content = content.unwrap();
        if let Some(tab_pos) = content.chars().position(|b| b == '\t') {
            let id = &content[..tab_pos];
            let id: u32 = id.parse().unwrap();
            writer
                .write_all(&id.to_le_bytes())
                .expect("unable to write id");

            total_records += 1;
        }
    }

    // Scratch buffer
    let mut seq_iupac: Vec<Iupac> = Vec::with_capacity(seq_len);

    // Read all sequences
    let reader = BufReader::new(File::open(input).expect("unable to open file"));
    for content in reader.lines() {
        let content = content.unwrap();
        if let Some(tab_pos) = content.chars().position(|b| b == '\t') {
            let seq = &content[tab_pos + 1..];
            for c in seq.chars() {
                seq_iupac.push(Iupac::from_utf8(c));
            }

            // SAFETY: Iupac is repr(u8)
            let bytes: &[u8] =
                unsafe { std::slice::from_raw_parts(seq_iupac.as_ptr().cast::<u8>(), seq_len) };
            writer.write_all(bytes).expect("unable to write sequence");

            seq_iupac.clear();
        }
    }

    let rename = format!("{}_{total_records}.bin", output.clone().to_string_lossy());
    std::fs::rename(output, rename).expect("unable to rename output file");
}
