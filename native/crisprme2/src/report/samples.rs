//! Sample presence sets for the assembly report.
//!
//! A row's sample payload is a set of per-sample haplotype-presence bitmasks,
//! interned to a `SampleSetId`. Two sets combine by *per-sample copy-wise OR*
//! (merge-join on sample), not plain set union, so `HG1:1|0` and `HG1:0|1`
//! collapse to `HG1:1|1` (homozygous). Ploidy is data-defined per sample; a
//! `u16` mask supports up to 16 copies, beyond any real assembly.

use ahash::AHashMap;

/// Per-sample copy-presence bitmask. Bit `b` set == present on copy `b`, where
/// copy `b` is `SampleTable::hap_layout[sample][b]` (declared-layout order, not
/// raw PanSN hap_id).
pub type Presence = u16;

/// Dense, ordered sample vocabulary plus per-sample copy layout. Built once
/// from `AssemblyInputs`; decodes ids and renders presence masks.
pub struct SampleTable {
    /// index -> sample name, sorted (stable u32 ids).
    names: Vec<String>,
    /// index -> that sample's PanSN hap_ids, sorted ascending. Bit position `b`
    /// of a presence mask maps to `hap_layout[sample][b]`. Length == ploidy.
    hap_layout: Vec<Vec<u32>>,
}

impl SampleTable {
    pub fn new(names: Vec<String>, hap_layout: Vec<Vec<u32>>) -> Self {
        debug_assert_eq!(names.len(), hap_layout.len());
        debug_assert!(hap_layout.iter().all(|h| h.len() <= Presence::BITS as usize));
        Self { names, hap_layout }
    }

    #[inline]
    pub fn name(&self, sample: u32) -> &str {
        &self.names[sample as usize]
    }

    #[inline]
    pub fn ploidy(&self, sample: u32) -> usize {
        self.hap_layout[sample as usize].len()
    }

    /// Bit position (0..ploidy) for a PanSN hap_id within a sample, or `None`
    /// if the hap_id is not part of this sample's declared layout.
    #[inline]
    pub fn bit_of(&self, sample: u32, hap_id: u32) -> Option<usize> {
        self.hap_layout[sample as usize].iter().position(|&h| h == hap_id)
    }

    /// Append one sample's presence mask as `1|0|1` over its declared copies.
    pub fn render_presence(&self, sample: u32, mask: Presence, out: &mut String) {
        for b in 0..self.ploidy(sample) {
            if b > 0 {
                out.push('|');
            }
            out.push(if mask & (1 << b) != 0 { '1' } else { '0' });
        }
    }
}

/// One canonical entry: a sample and its copy-presence bitmask.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct SamplePresence {
    pub sample: u32,
    pub mask: Presence,
}

/// Canonical sample set: at most one entry per sample, sorted by sample id.
/// Canonicalization (fold duplicate samples by OR, sort) happens on intern.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct SampleSet(Box<[(u32, Presence)]>);

impl SampleSet {
    /// The canonical `(sample, mask)` entries, sorted by sample id.
    #[inline]
    pub fn entries(&self) -> &[(u32, Presence)] {
        &self.0
    }
}

/// Interned handle into `SampleSetRegistry`. `u32` indexes distinct sets — not
/// a 32-sample cap.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct SampleSetId(pub u32);

pub struct SampleSetRegistry {
    sets: Vec<SampleSet>,
    intern: AHashMap<SampleSet, SampleSetId>,
}

impl Default for SampleSetRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl SampleSetRegistry {
    pub fn new() -> Self {
        Self { sets: Vec::new(), intern: AHashMap::new() }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.sets.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sets.is_empty()
    }

    /// Intern raw `(sample, mask)` pairs. Folds duplicate samples by OR and
    /// sorts, so any permutation or pre-split of the same biological content
    /// interns to the same id (`HG1:1|0` + `HG1:0|1` == `HG1:1|1`). This is the
    /// single place same-sample copies are folded.
    pub fn intern(&mut self, pairs: &[SamplePresence]) -> SampleSetId {
        let mut folded: AHashMap<u32, Presence> = AHashMap::with_capacity(pairs.len());
        for p in pairs {
            *folded.entry(p.sample).or_insert(0) |= p.mask;
        }
        let mut v: Vec<(u32, Presence)> = folded.into_iter().collect();
        v.sort_unstable_by_key(|&(s, _)| s);
        let set = SampleSet(v.into_boxed_slice());
        if let Some(&id) = self.intern.get(&set) {
            return id;
        }
        let id = SampleSetId(self.sets.len() as u32);
        self.sets.push(set.clone());
        self.intern.insert(set, id);
        id
    }

    /// Union two sets by per-sample copy-wise OR (merge-join on sample), then
    /// intern. The D3 cross-assembly merge primitive; the OR-collapse of a
    /// shared sample is automatic via `intern`'s fold.
    pub fn union(&mut self, a: SampleSetId, b: SampleSetId) -> SampleSetId {
        let sa = &self.sets[a.0 as usize].0;
        let sb = &self.sets[b.0 as usize].0;
        let mut merged: Vec<SamplePresence> = Vec::with_capacity(sa.len() + sb.len());
        merged.extend(sa.iter().map(|&(s, m)| SamplePresence { sample: s, mask: m }));
        merged.extend(sb.iter().map(|&(s, m)| SamplePresence { sample: s, mask: m }));
        self.intern(&merged)
    }

    /// Append `HG1:1|1,HG2:1|0` for the report `samples` column — samples in id
    /// order, copies in declared-layout order.
    pub fn render(&self, id: SampleSetId, table: &SampleTable, out: &mut String) {
        for (i, &(sample, mask)) in self.sets[id.0 as usize].0.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(table.name(sample));
            out.push(':');
            table.render_presence(sample, mask, out);
        }
    }

    /// Canonical `(sample, mask)` entries of an interned set, sorted by sample
    #[inline]
    pub fn entries_of(&self, id: SampleSetId) -> &[(u32, Presence)] {
        &self.sets[id.0 as usize].0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table() -> SampleTable {
        // HG1 diploid (haps 1,2), HG2 tetraploid (haps 0,1,2,3)
        SampleTable::new(
            vec!["HG1".into(), "HG2".into()],
            vec![vec![1, 2], vec![0, 1, 2, 3]],
        )
    }
    fn sp(sample: u32, mask: Presence) -> SamplePresence {
        SamplePresence { sample, mask }
    }

    #[test]
    fn homozygous_fold_on_intern() {
        let (t, mut r) = (table(), SampleSetRegistry::new());
        // HG1 present on copy 0 and (separately) copy 1 -> folds to 1|1
        let a = r.intern(&[sp(0, 0b01), sp(0, 0b10)]);
        let b = r.intern(&[sp(0, 0b11)]);
        assert_eq!(a, b); // same interned id
        let mut s = String::new();
        r.render(a, &t, &mut s);
        assert_eq!(s, "HG1:1|1");
    }

    #[test]
    fn union_is_per_sample_or() {
        let (t, mut r) = (table(), SampleSetRegistry::new());
        let a = r.intern(&[sp(0, 0b10), sp(1, 0b0001)]); // HG1:0|1, HG2:1|0|0|0
        let b = r.intern(&[sp(0, 0b01)]);                // HG1:1|0
        let u = r.union(a, b);
        let mut s = String::new();
        r.render(u, &t, &mut s);
        assert_eq!(s, "HG1:1|1,HG2:1|0|0|0"); // HG1 OR-collapsed; HG2 carried
    }

    #[test]
    fn render_uses_declared_ploidy() {
        let (t, mut r) = (table(), SampleSetRegistry::new());
        let id = r.intern(&[sp(1, 0b1010)]); // HG2 copies 1 and 3
        let mut s = String::new();
        r.render(id, &t, &mut s);
        assert_eq!(s, "HG2:0|1|0|1");
    }

    #[test]
    fn bit_of_maps_pansn_hapid_to_position() {
        let t = table();
        assert_eq!(t.bit_of(0, 1), Some(0)); // HG1 hap_id 1 -> bit 0
        assert_eq!(t.bit_of(0, 2), Some(1)); // HG1 hap_id 2 -> bit 1
        assert_eq!(t.bit_of(0, 0), None);    // hap_id 0 not in HG1's layout
    }
}