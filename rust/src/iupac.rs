pub fn matches_iupac(nt: u8, pattern: u8) -> bool {
    const MAP: &[(u8, &[u8])] = &[
        (b'A', &[b'A']), (b'C', &[b'C']), (b'G', &[b'G']), (b'T', &[b'T']),
        (b'R', &[b'A', b'G']), (b'Y', &[b'C', b'T']),
        (b'S', &[b'G', b'C']), (b'W', &[b'A', b'T']),
        (b'K', &[b'G', b'T']), (b'M', &[b'A', b'C']),
        (b'B', &[b'C', b'G', b'T']), (b'D', &[b'A', b'G', b'T']),
        (b'H', &[b'A', b'C', b'T']), (b'V', &[b'A', b'C', b'G']),
        (b'N', &[b'A', b'C', b'G', b'T']),
    ];

    MAP.iter()
        .find(|(key, _)| *key == pattern)
        .map(|(_, allowed)| allowed.contains(&nt))
        .unwrap_or(false)
}