/*
#[test]
fn correct_results_against_bio_crate() {
    const FASTA_SIZE: usize = 100000;
    const SEQ_LEN: usize = 24;
    const GUIDE_LEN: usize = 20;
    const DELTA: usize = 1;
    const THRESHOLD: u8 = 6;

    // For debug
    let score_fn = |a: u8, b: u8| if a == b { 0i32 } else { -1i32 };
    let scoring = Scoring::new(0, -1, &score_fn).xclip(MIN_SCORE).yclip(0);

    // Semi-global aligner
    let mut aligner = Aligner::with_scoring(scoring);

    // Load reference fasta file
    let fasta = load_from_file("../fasta/chr22.fa").expect("fasta file not found");

    // Test a subset of the fasta file
    let fasta = fasta[0..FASTA_SIZE].to_vec();
    let (windows, n) = materialize_windows(&fasta, SEQ_LEN, DELTA);
    let guide = fasta[DELTA * 20..DELTA * 20 + GUIDE_LEN].to_vec();

    // Run CPU aligner
    let mut cpu_valid = 0;
    for i in 0..n {
        let beg = i * DELTA;
        let aligment = aligner.custom(&guide, &fasta[beg..beg + SEQ_LEN]);
        if (-aligment.score) as u8 <= THRESHOLD {
            cpu_valid += 1;
        }
    }

    let result = cuda::filter(&guide, &windows, SEQ_LEN, n);
    let gpu_valid = result
        .iter()
        .map(|e| if *e <= THRESHOLD { 1 } else { 0 })
        .sum::<u32>();

    // We must find the same number of valid sequences
    assert_eq!(cpu_valid, gpu_valid);

    // Store all below-threshold sequences
    let result: Vec<(Vec<u8>, u8)> = result
        .iter()
        .enumerate()
        .filter_map(|(i, score)| {
            if *score <= THRESHOLD {
                let beg = i * DELTA;
                let seq = fasta[beg..beg + SEQ_LEN].to_vec();
                Some((seq, *score))
            } else {
                None
            }
        })
        .collect();

    // Check that the GPU result is identical to the CPU result
    for (seq, gpu_score) in &result {
        let aligment = aligner.custom(&guide, seq);
        assert_eq!(*gpu_score, (-aligment.score) as u8);
    }
}
*/
