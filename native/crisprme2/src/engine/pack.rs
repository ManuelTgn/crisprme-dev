use crate::batching::batching::{unpack_occ, WindowBatch};
use crate::memory::batch::WindowRingBatch;
use crate::sequence::iupac::Iupac;

/*
fn fill_window_ring(dst: &mut SequenceRingBatch, src: WindowBatch, sequence_len: usize) {
    debug_assert_eq!(sequence_len, dst.descriptor.sequence_len);
    debug_assert_eq!(src.windows.len(), src.occs.len());

    // Set actual counts for this slot
    dst.descriptor.window_count = src.windows.len();
    dst.descriptor.occ_count = src.total_hits;

    // --- 1) Copy unique windows (packed at start of slot) ---
    let slen = sequence_len;
    let out = dst.windows_iupac_mut(); // len = window_count * slen

    for (i, w) in src.windows.iter().enumerate() {
        let base = i * slen;
        for j in 0..slen {
            out[base + j] = Iupac::new(w[j]);
        }
    }

    // --- 2) Flatten occurrences into parallel arrays + prefix index ---
    let starts = dst.occ_starts_mut();
    let lens   = dst.occ_lens_mut();
    let contig = dst.occ_contig_mut();
    let pos    = dst.occ_pos_mut();
    let strand = dst.occ_strand_mut();

    let mut cursor: u32 = 0;
    for (i, occs) in src.occs.iter().enumerate() {
        starts[i] = cursor;
        lens[i] = occs.len() as u32;

        for &occ in occs {
            let (cid, p, s) = unpack_occ(occ);
            let k = cursor as usize;
            contig[k] = cid;
            pos[k] = p;
            strand[k] = s; // 1 = '+', 0 = '-'
            cursor += 1;
        }
    }

    // --- 3) Only windows need to go to GPU for mining ---
    dst.sync_windows_cpu_to_gpu(); // sync only window bytes
}
    */
