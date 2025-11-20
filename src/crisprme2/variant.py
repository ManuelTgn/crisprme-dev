""" """

from .crisprme2_error import Crisprme2VariantRecordError
from .logger import CrisprmeLoggers
from .sample import Sample

from typing import List, Set

from cyvcf2 import Variant
from time import time

import numpy as np

import os

# variants types
VTYPES = ["snv", "indel"]

class VariantRecord:

    def __init__(self, variant: Variant, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._variant = variant  # read vcf line 
        self._allelesnum = len(self.alt)
        self._vtype = self._assign_variant_type()  # compute variant types
        self._vid = self._assign_id()  # compute variant ids        


    def __repr__(self) -> str:
        altalleles = ",".join(self.alt)
        return (
            f'<{self.__class__.__name__} object; variant="{self.contig} {self.position} '
            f'{self.ref} {altalleles}">'
        )

    def __str__(self) -> str:
        altalleles = ",".join(self.alt)
        return f"{self.contig}\t{self.position}\t{self.ref}\t{altalleles}"

    def __eq__(self, vrecord: object) -> bool:
        if not isinstance(vrecord, VariantRecord):
            return NotImplemented
        if not hasattr(vrecord, "_variant"):  
            self._loggers.errorlog.log_raise_exception(f"Comparison between {self.__class__.__name__} object failed", os.EX_DATAERR, Crisprme2VariantRecordError)
        return (
            self._variant.CHROM == vrecord.contig
            and self._variant.POS == vrecord.position
            and self._variant.REF == vrecord.ref
            and self._variant.ALT == vrecord.alt
        )

    def __lt__(self, vrecord: "VariantRecord") -> bool:
        return self._variant.POS < vrecord.position

    def __gt__(self, vrecord: "VariantRecord") -> bool:
        """Compare two VariantRecord objects based on position.

        Compares the current VariantRecord object with another VariantRecord object
        based on their positions.

        Args:
            vrecord: The other VariantRecord object to compare to.

        Returns:
            True if the current object's position is greater than the other object's
            position, False otherwise.
        """
        return self._variant.POS > vrecord.position

    def __hash__(self) -> int:
        """Return a hash value for the variant record.

        Returns a hash value for the VariantRecord object, based on its chromosome,
        position, reference allele, and alternative alleles. This allows
        VariantRecord objects to be used as keys in dictionaries or sets.

        Returns:
            The hash value of the variant record.
        """
        return hash((self._variant.CHROM, self._variant.POS, self._variant.REF, tuple(self._variant.ALT)))


    def _assign_id(self) -> List[str]:
        if self._allelesnum == 1:
            # variant id not available, construct the id using chrom, position, ref,
            # and alt (e.g. chrx-100-A/G)
            return [_compute_id(self._variant.CHROM, self._variant.POS, self._variant.REF, self._variant.ALT[0])]
        # if multiallelic site compute the id for each alternative allele
        # avoid potential confusion due to alternative alleles at same position
        # labeled with the same id
        return [
            _compute_id(self._variant.CHROM, self._variant.POS, self._variant.REF, altallele)
            for altallele in self._variant.ALT
        ]
    
    def _assign_variant_type(self) -> List[str]:
        return [VTYPES[0] if len(self.ref) == len(altallele) else VTYPES[1] for altallele in self.alt]

            
    @property
    def contig(self) -> str:
        return self._variant.CHROM
    
    @property
    def position(self) -> int:
        return self._variant.POS
    
    @property
    def ref(self) -> str:
        return self._variant.REF

    @property
    def alt(self) -> List[str]:
        return self._variant.ALT
    
    @property
    def filters(self) -> List[str]:
        return self._variant.FILTERS
    
    @property
    def id(self) -> List[str]:
        return self._vid
    
    @property
    def af(self) -> List[float]:
        return self._variant.aaf

    @property
    def ploidy(self) -> int:
        return self._variant.ploidy

    @property
    def genotypes(self) -> List[List[int]]:
        return self._variant.genotypes

    @property
    def num_called(self) -> int:
        return self._variant.num_called

    @property
    def num_hom_ref(self) -> int:
        return self._variant.num_hom_ref
    
    @property
    def vtype(self) -> List[str]:
        return self._vtype
    


#     def _copy(self, i: int) -> "VariantRecord":
#         """Create a copy of the variant record for a specific allele.

#         Creates a copy of the current VariantRecord object, representing a single
#         alternative allele specified by the index `i`.  This is useful for
#         handling multiallelic sites where each alternative allele needs to be
#         treated as a separate variant.

#         Args:
#             i: The index of the alternative allele to copy.

#         Returns:
#             A new VariantRecord object representing the specified alternative allele.
#         """
#         # sourcery skip: class-extract-method
#         # copy current variant record instance
#         vrecord = VariantRecord(self._debug)  # create new instance
#         # adjust ref/alt alleles and positions for multiallelic sites
#         ref, alt, position = adjust_multiallelic(
#             self._ref, self._alt[i], self._position
#         )
#         vrecord._chrom = self._chrom
#         vrecord._position = position
#         vrecord._ref = ref
#         vrecord._alt = [alt]
#         vrecord._allelesnum = 1
#         vrecord._vtype = [self._vtype[i]]
#         vrecord._filter = self._filter
#         vrecord._afs = [self._afs[i]]
#         vrecord._vid = [self._vid[i]]
#         vrecord._samples = [self._samples[i]]
#         return vrecord

#     def split(self, vtype: Optional[str] = None) -> List["VariantRecord"]:
#         """Split a multiallelic variant record by variant type.

#         Splits a multiallelic VariantRecord object into a list of VariantRecord objects,
#         each representing a single alternative allele of the specified variant type.

#         Args:
#             vtype: The variant type to select ("snp" or "indel"). If None, all
#                 variant types are included.

#         Returns:
#             A list of VariantRecord objects, one for each alternative allele matching
#             the specified variant type.
#         """
#         vtypes_filter = VTYPES if vtype is None else [vtype]
#         return [
#             self._copy(i)
#             for i, _ in enumerate(self._vtype)
#             if self._vtype[i] in vtypes_filter
#         ]

#     def get_altalleles(self, vtype: str) -> List[str]:
#         """Retrieve alternative alleles of a specific variant type.

#         Returns a list of alternative alleles that match the specified variant type
#         (either "snp" or "indel").

#         Args:
#             vtype: The variant type to select ("snp" or "indel").

#         Returns:
#             A list of alternative alleles matching the specified type.
#         """
#         assert vtype in VTYPES
#         # return the alternative alleles representing snps or indels
#         return [
#             altallele
#             for i, altallele in enumerate(self._alt)
#             if self._vtype[i] == vtype
#         ]

#     def pytest_initialize(
#         self, position: int, ref: str, alt: str, vtype: str, vid: str, afs: List[float]
#     ) -> None:
#         """Initialize the VariantRecord for pytest with provided values.

#         Sets the chromosome, position, reference allele, alternative allele,
#         variant type, variant ID, and allele frequencies for testing purposes.

#         Args:
#             position (int): The variant position.
#             ref (str): The reference allele.
#             alt (str): The alternative allele.
#             vtype (str): The variant type.
#             vid (str): The variant ID.
#             afs (List[float]): The allele frequencies.
#         """
#         self._chrom = "chrx"
#         self._position = position
#         self._ref = ref
#         self._alt = [alt]
#         self._vtype = [vtype]
#         self._vid = [vid]
#         self._afs = afs

#     @property
#     def filter(self) -> str:
#         return self._filter

#     @property
#     def contig(self) -> str:
#         return self._chrom

#     @property
#     def position(self) -> int:
#         return self._position

#     @property
#     def ref(self) -> str:
#         return self._ref

#     @property
#     def alt(self) -> List[str]:
#         return self._alt

#     @property
#     def vtype(self) -> List[str]:
#         return self._vtype

#     @property
#     def afs(self) -> List[float]:
#         return self._afs

#     @property
#     def samples(self) -> List[Tuple[Set[str], Set[str]]]:
#         return self._samples

#     @property
#     def id(self) -> List[str]:
#         return self._vid

#     @property
#     def allelesnum(self) -> int:
#         return self._allelesnum


def _assign_vtype(ref: str, alt: str) -> str:
    return VTYPES[1] if len(ref) != len(alt) else VTYPES[0]


def _compute_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    # compute variant id for variants without id, or multiallelic sites
    # use IGVF consortium notation
    return f"{chrom}-{pos}-{ref}/{alt}"


# def adjust_multiallelic(ref: str, alt: str, pos: int) -> Tuple[str, str, int]:
#     """Adjust reference/alternative alleles and position for multiallelic sites.

#     Adjusts the reference and alternative alleles, and the variant position for
#     multiallelic sites based on the lengths of the original reference and
#     alternative alleles.  This function helps normalize variant representation
#     for easier comparison and processing. The function assumes multiallelic
#     variants are left-aligned.

#     Args:
#         ref: The original reference allele.
#         alt: The original alternative allele.
#         pos: The original variant position.

#     Returns:
#         A tuple containing the adjusted reference allele, alternative allele, and
#         variant position.
#     """

#     if len(ref) == len(alt):  # likely snp
#         ref_new, alt_new = ref[0], alt[0]  # adjust ref/alt alleles
#         pos_new = pos  # ref/alt have same length
#     elif len(ref) > len(alt):  # deletion
#         ref_new = ref[len(alt) - 1 :]  # adjust ref allele
#         alt_new = alt[-1]  # adjust alt allele
#         pos_new = pos + (len(alt)) - 1  # adjust variant position
#     else:  # insertion
#         ref_new = ref[-1]  # adjust ref allele
#         alt_new = alt[len(ref) - 1 :]  # adjust alt allele
#         pos_new = pos + len(ref) - 1  # adjust variant position
#     return ref_new, alt_new, pos_new





# def _parse_genotype_unphased(
#     gt_alleles: List[str],
#     sample: str,
#     sampleshap: List[Tuple[Set[str], Set[str]]],
# ) -> List[Tuple[Set[str], Set[str]]]:
#     """Parse unphased genotype alleles and update sample sets.

#     Updates the sample sets for each alternative allele based on the unphased
#     genotype information. Assigns the sample to the appropriate set for each
#     allele present in the genotype.

#     Args:
#         gt_alleles: A list containing the alleles for the two haplotypes or more.
#         sample: The sample name.
#         sampleshap: A list of tuples of sets, tracking samples for each allele
#             and haplotype.

#     Returns:
#         The updated list of tuples of sets with sample assignments.
#     """
#     if len(gt_alleles) != 2:  # handle genotypes like 0/1/2
#         for gt in gt_alleles:
#             if gt not in ["0", "."]:
#                 sampleshap[int(gt) - 1][0].add(sample)
#     else:  # handle genotypes like 0/1
#         gt1, gt2 = gt_alleles  # retrieve allele occurring on first and second copy
#         if gt1 not in ["0", "."] and gt1 == gt2:  # special case 1/1
#             sampleshap[int(gt1) - 1][0].add(sample)
#             sampleshap[int(gt2) - 1][1].add(sample)
#         else:
#             if gt1 not in ["0", "."]:  # 1/0
#                 sampleshap[int(gt1) - 1][0].add(sample)
#             if gt2 not in ["0", "."]:  # 0/1
#                 sampleshap[int(gt2) - 1][0].add(sample)
#     return sampleshap



def _split_vcfline(vcfline: str, loggers: CrisprmeLoggers) -> List[str]:
    try:
        return vcfline.strip().split()
    except Exception as e:
        loggers.errorlog.log_exception(f"VCF line parsing failed: {e}", os.EX_DATAERR)
