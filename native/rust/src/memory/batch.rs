use std::ops::Deref;

use bincode::de;

use super::ring::{RingAdapter, RingSlotLease};

// use crate::bindings;
use crate::crispr::guide::Guide;
use crate::sequence::iupac::Iupac;
use crate::sequence::sequence::Sequence;
use crate::alignment::alignment::Alignment;
use crate::storage::reader::SequenceBatchDescr;
use crate::storage::writer::AlignmentBatchDescr;

/// View into a batch of sequences in a ring buffer
pub struct SequenceRingBatch {
    pub descriptor: SequenceBatchDescr,
    pub lease: RingSlotLease,
}

impl Deref for SequenceRingBatch {
    type Target = RingSlotLease;
    fn deref(&self) -> &Self::Target {
        &self.lease
    }
}

/// Allow to attach to ring buffer
impl RingAdapter for SequenceRingBatch {
    type Descr = SequenceBatchDescr;

    fn attach(lease: RingSlotLease, descr: Self::Descr) -> Self {
        Self {
            descriptor: descr,
            lease,
        }
    }

    fn detach(self) -> (Self::Descr, RingSlotLease) {
        (self.descriptor, self.lease)
    }

    fn as_mut(&mut self) -> &mut RingSlotLease {
        &mut self.lease
    }

    fn as_ref(&self) -> &RingSlotLease {
        &self.lease
    }
}

impl SequenceRingBatch {

    /// Apply a mask to filter sequences, this compacts only valid sequences to the start of the batch
    /// and changes the size of the batch. No memory is deallocated in this operation.
    pub fn apply_mask(&mut self, mask: &[bool], sync: bool) {
        assert_eq!(
            self.descriptor.sequence_count,
            mask.len(),
            "mask must have the same lenght as the batch!"
        );

        // Copy memory from GPU if necessary
        if sync {
            self.lease.sync_gpu_to_cpu(None);
        }

        let slen = self.descriptor.sequence_len;

        // Mask the sequences
        let mut w = 0;
        let sequences = self.iupac_mut();
        for (r, keep) in mask.iter().enumerate() {
            if *keep {
                if r != w {
                    sequences.copy_within(
                        r * slen..(r + 1) * slen,
                        w * slen,
                    );
                }
                w += 1;
            }
        }

        // Mask the ids
        w = 0;
        let ids = &mut self.ids_mut();
        for (r, keep) in mask.iter().enumerate() {
            if *keep {
                if r != w {
                    ids[w] = ids[r];
                }
                w += 1;
            }
        }

        // Update batch size
        self.descriptor.sequence_count = w;

        // Sync GPU memory
        if sync {
            self.lease.sync_cpu_to_gpu(None);
        }
    }

    /// Returns number of sequences
    pub fn len(&self) -> usize {
        self.descriptor.sequence_count
    }

    /// Returns pointer to GPU memory
    pub fn gpu_ptr(&self) -> *const u8 {
        self.lease.gpu_ptr()
    }

    // /// Calculate the edit-distance score between all sequences and the guide
    // pub fn edit_distace_scores(&mut self, guide: &Guide, result: &mut [u8]) {
    //     // TODO: This allocates on gpu every time
    //     bindings::scores::scores_into(
    //         guide.as_slice(),
    //         self.iupac(),
    //         result,
    //         self.descriptor.sequence_len,
    //         self.descriptor.sequence_count,
    //     );
    // }

    /// Iterator over all sequences
    pub fn sequences(&self) -> impl Iterator<Item = Sequence<'_>> {
        self.iupac()
            .chunks_exact(self.descriptor.sequence_len)
            .map(Sequence::new)
    }

    /// Iterator over all sequences and ids
    pub fn sequences_with_ids(&self) -> impl Iterator<Item = (u32, Sequence<'_>)> {
        self.ids().iter().cloned().zip(self.sequences())
    }

    /// Mutable slice of iupac data
    pub fn iupac_mut(&mut self) -> &mut [Iupac] {
        unsafe {
            // SAFETY: Iupac is repr(u8)
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr() as *mut Iupac,
                self.descriptor.sequence_count * self.descriptor.sequence_len,
            )
        }
    }

    /// Slice of iupac data
    pub fn iupac(&self) -> &[Iupac] {
        unsafe {
            // SAFETY: Iupac is repr(u8)
            std::slice::from_raw_parts(
                self.lease.as_ptr() as *const Iupac,
                self.descriptor.sequence_count * self.descriptor.sequence_len,
            )
        }
    }

    /// Mutable slice of sequence ids
    pub fn ids_mut(&mut self) -> &mut [u32] {
        let alignment = align_of::<u32>();
        let offset = (self.iupac().len() + alignment - 1) & !(alignment - 1);
        unsafe {
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr().add(offset) as *mut u32,
                self.descriptor.sequence_count,
            )
        }
    }

    /// Slice of sequence ids
    pub fn ids(&self) -> &[u32] {
        let alignment = align_of::<u32>();
        let offset = (self.iupac().len() + alignment - 1) & !(alignment - 1);
        unsafe {
            std::slice::from_raw_parts(
                self.lease.as_ptr().add(offset) as *const u32,
                self.descriptor.sequence_count,
            )
        }
    }
}

/// View into a batch of alignments inside a ring buffer
pub struct AlignmentRingBatch {
    descriptor: AlignmentBatchDescr,
    lease: RingSlotLease,
}

impl AlignmentRingBatch {

    /// Return pointer to GPU memory
    pub fn gpu_ptr_mut(&mut self) -> *mut u8 {
        self.lease.gpu_ptr_mut()
    }

    /// Returns the total available alignments
    pub fn capacity(&self) -> usize {
        self.lease.capacity() / size_of::<Alignment>()
    }

    /// Returns number of alignments
    pub fn len(&self) -> usize {
        self.descriptor.alignment_count
    }

    /// Set the amount of mined alignments
    pub fn set_len(&mut self, len: usize) {
        assert!(len < self.capacity());
        self.descriptor.alignment_count = len;
    }

    pub fn alignments(&self) -> &[Alignment] {
        unsafe {
            std::slice::from_raw_parts(
                self.lease.as_ptr() as *const Alignment,
                self.descriptor.alignment_count,
            )
        }
    }

    pub fn alignments_mut(&mut self) -> &mut [Alignment] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr() as *mut Alignment,
                self.descriptor.alignment_count,
            )
        }
    }

    pub fn replace_pos_by_id(&mut self, batch: &SequenceRingBatch) {
        let ids = batch.ids();
        for align in self.alignments_mut() {
            align.id = ids[align.id as usize]; 
        }
    }

    #[inline(always)]
    pub fn sync_gpu_to_cpu(&mut self, bytes: Option<usize>) {
        self.lease.sync_gpu_to_cpu(bytes);
    }

    #[inline(always)]
    pub fn sync_cpu_to_gpu(&mut self, bytes: Option<usize>) {
        self.lease.sync_cpu_to_gpu(bytes);
    }
}

impl RingAdapter for AlignmentRingBatch {
    type Descr = AlignmentBatchDescr;

    fn attach(lease: RingSlotLease, descr: Self::Descr) -> Self {
        Self {
            descriptor: descr,
            lease,
        }
    }

    fn detach(self) -> (Self::Descr, RingSlotLease) {
        (self.descriptor, self.lease)
    }

    fn as_mut(&mut self) -> &mut RingSlotLease {
        &mut self.lease
    }

    fn as_ref(&self) -> &RingSlotLease {
        &self.lease
    }
}




#[derive(Clone, Copy, Default, Debug)]
pub struct WindowBatchDescr {
    // size (spacer+pam window length)
    pub sequence_len: usize,  
    // number of unique windows
    pub window_count: usize,  
    // total occurrences across all windows
    pub occ_count: usize,
}

pub struct WindowRingBatch {
    pub descriptor: WindowBatchDescr,
    pub lease: RingSlotLease,
}

impl Deref for WindowRingBatch {
    type Target = RingSlotLease;
    fn deref(&self) -> &Self::Target { &self.lease }
}

impl RingAdapter for WindowRingBatch {
    type Descr = WindowBatchDescr;

    fn attach(lease: RingSlotLease, descr: Self::Descr) -> Self {
        Self { descriptor: descr, lease }
    }

    fn detach(self) -> (Self::Descr, RingSlotLease) {
        (self.descriptor, self.lease)
    }

    fn as_mut(&mut self) -> &mut RingSlotLease { &mut self.lease }

    fn as_ref(&self) -> &RingSlotLease { &self.lease }
}

impl WindowRingBatch {
    #[inline] pub fn windows_len(&self) -> usize { self.descriptor.window_count }

    #[inline] pub fn occ_len(&self) -> usize { self.descriptor.occ_count }

    #[inline] fn windows_bytes_len(&self) -> usize {
        self.descriptor.window_count * self.descriptor.sequence_len
    }

    #[inline] fn off_windows_end(&self) -> usize {
        self.windows_bytes_len()
    }

    #[inline] fn off_occ_starts(&self) -> usize {
        // align to u32
        let a = align_of::<u32>();
        (self.off_windows_end() + a - 1) & !(a - 1)
    }

    #[inline] fn off_occ_lens(&self) -> usize {
        self.off_occ_starts() + self.descriptor.window_count * size_of::<u32>()
    }

    #[inline] fn off_occ_contig(&self) -> usize {
        self.off_occ_lens() + self.descriptor.window_count * size_of::<u32>()
    }

    #[inline] fn off_occ_pos(&self) -> usize {
        self.off_occ_contig() + self.descriptor.occ_count * size_of::<u32>()
    }

    #[inline] fn off_occ_strand(&self) -> usize {
        self.off_occ_pos() + self.descriptor.occ_count * size_of::<u32>()
    }

    /// The bytes that must be synced to GPU for mining (windows only).
    #[inline] pub fn gpu_windows_bytes(&self) -> usize {
        // If your CUDA miner expects sequences at start of slot: perfect.
        self.windows_bytes_len()
    }

    /// Windows (unique) as IUPAC slice
    pub fn windows_iupac(&self) -> &[Iupac] {
        unsafe {
            std::slice::from_raw_parts(
                self.lease.as_ptr() as *const Iupac,
                self.windows_bytes_len(),
            )
        }
    }

    pub fn windows_iupac_mut(&mut self) -> &mut [Iupac] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr() as *mut Iupac,
                self.windows_bytes_len(),
            )
        }
    }

    pub fn occ_starts(&self) -> &[u32] {
        unsafe {
            std::slice::from_raw_parts(
                self.lease.as_ptr().add(self.off_occ_starts()) as *const u32,
                self.descriptor.window_count,
            )
        }
    }

    pub fn occ_starts_mut(&mut self) -> &mut [u32] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr().add(self.off_occ_starts()) as *mut u32,
                self.descriptor.window_count,
            )
        }
    }

    pub fn occ_lens(&self) -> &[u32] {
        unsafe {
            std::slice::from_raw_parts(
                self.lease.as_ptr().add(self.off_occ_lens()) as *const u32,
                self.descriptor.window_count,
            )
        }
    }

    pub fn occ_lens_mut(&mut self) -> &mut [u32] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr().add(self.off_occ_lens()) as *mut u32,
                self.descriptor.window_count,
            )
        }
    }

    pub fn occ_contig(&self) -> &[u32] {
        unsafe {
            std::slice::from_raw_parts(
                self.lease.as_ptr().add(self.off_occ_contig()) as *const u32,
                self.descriptor.occ_count,
            )
        }
    }

    pub fn occ_contig_mut(&mut self) -> &mut [u32] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr().add(self.off_occ_contig()) as *mut u32,
                self.descriptor.occ_count,
            )
        }
    }

    pub fn occ_pos(&self) -> &[u32] {
        unsafe {
            std::slice::from_raw_parts(
                self.lease.as_ptr().add(self.off_occ_pos()) as *const u32,
                self.descriptor.occ_count,
            )
        }
    }

    pub fn occ_pos_mut(&mut self) -> &mut [u32] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr().add(self.off_occ_pos()) as *mut u32,
                self.descriptor.occ_count,
            )
        }
    }

    pub fn occ_strand(&self) -> &[u8] {
        unsafe {
            std::slice::from_raw_parts(
                self.lease.as_ptr().add(self.off_occ_strand()) as *const u8,
                self.descriptor.occ_count,
            )
        }
    }

    pub fn occ_strand_mut(&mut self) -> &mut [u8] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.lease.as_mut_ptr().add(self.off_occ_strand()) as *mut u8,
                self.descriptor.occ_count,
            )
        }
    }

    #[inline(always)]
    pub fn sync_windows_cpu_to_gpu(&mut self) {
        let bytes = self.gpu_windows_bytes();
        self.lease.sync_cpu_to_gpu(Some(bytes));
    }
}

pub fn window_ring_slot_bytes(sequence_len: usize, max_windows: usize, max_occs: usize) -> usize {
    let windows = max_windows * sequence_len;               // u8
    let pad = align_of::<u32>();                            // padding
    let starts = max_windows * size_of::<u32>();
    let lens   = max_windows * size_of::<u32>();
    let contig = max_occs   * size_of::<u32>();
    let pos    = max_occs   * size_of::<u32>();
    let strand = max_occs   * size_of::<u8>();
    windows + pad + starts + lens + contig + pos + strand
}





// =================================================================================
// STD implementations

impl std::fmt::Debug for SequenceRingBatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "SequenceBatch(len: {}, data: [", self.len())?;
        for (id, seq) in self.sequences_with_ids() {
            writeln!(f, "\t(id: {id}, seq: {seq:?})")?;
        }
        writeln!(f, "])")?;
        Ok(())
    }
}

impl std::fmt::Debug for AlignmentRingBatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "AlignmentBatch(len: {}, data: [", self.len())?;
        for alig in self.alignments() {
            writeln!(f, "\t{alig:?}")?;
        }
        writeln!(f, "])")?;
        Ok(())
    }
}
