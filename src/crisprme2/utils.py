"""
Utility functions and constants for the CRISPRme2 tool.

This module provides helper functions for file and directory management, sequence
manipulation, IUPAC matching, and model extraction. It also defines shared constants
and static variables used across the CRISPRme2 software.
"""

from typing import List, Any
from itertools import permutations

import os

# ==============================================================================
#
# STATIC VARIABLES
#
# ==============================================================================

# dna alphabet
DNA = ["A", "C", "G", "T", "N"]

# complete iupac alphabet
IUPAC = DNA + ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]

# reverse complement dictionary
RC = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "M": "K",
    "K": "M",
    "H": "D",
    "D": "H",
    "B": "V",
    "V": "B",
    "N": "N",
    "S": "S",
    "W": "W",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "u": "a",
    "r": "y",
    "y": "r",
    "m": "k",
    "k": "m",
    "h": "d",
    "d": "h",
    "b": "v",
    "v": "b",
    "n": "n",
    "s": "s",
    "w": "w",
}

# dictionary to encode nucleotides combinations as iupac characters
IUPACTABLE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "H": "ACT",
    "B": "CGT",
    "V": "ACG",
    "D": "AGT",
    "N": "ACGT",
}

# dictionary to encode nucleotide strings as iupac characters
IUPAC_ENCODER = {
    perm: k
    for k, v in IUPACTABLE.items()
    for perm in {"".join(p) for p in permutations(v)}
}

# define strand directions: 0 -> 5'-3'; 1 -> 3'-5'
STRAND = [0, 1]

# tbi index file extension
TBI = "tbi"

# fai index file extension
FAI = "fai"

# off-targets length
OFFTARGETLEN = 30

# define VCF extensions
# vcf extensions
VCFEXTENSIONS = {"vcf", "vcf.gz", "bcf", "bcf.gz"}

# ==============================================================================
#
# UTILS MODULE FUNCTIONS
#
# ==============================================================================


def flatten_list(lst: List[List[Any]]) -> List[Any]:
    """Flatten a list of lists into a single list.

    Combines all elements from nested lists into a single flat list.

    Args:
        lst (List[List[Any]]): The list of lists to flatten.

    Returns:
        List[Any]: The flattened list.
    """
    return [e for sublist in lst for e in sublist]


def find_tbi_index(fname: str) -> bool:
    """Check if a Tabix index exists for the input VCF/BED file.

    Checks if a Tabix index (.tbi) exists for the given VCF/BED file and is a
    non-empty file.

    Args:
        fname: The path to the VCF/BED file.

    Returns:
        True if the index exists and is a non-empty file, False otherwise.
    """
    # avoid unexpected crashes due to file location
    tbi_index = f"{os.path.abspath(fname)}.{TBI}"
    if os.path.exists(tbi_index):  # index must be a non empty file
        return os.path.isfile(tbi_index) and os.stat(tbi_index).st_size > 0
    return False


def find_fai_index(fname: str) -> bool:
    # avoid unexpected crashes due to file location
    fai_index = f"{os.path.abspath(fname)}.{FAI}"
    if os.path.exists(fai_index):  # index must be a non empty file
        return os.path.isfile(fai_index) and os.stat(fai_index).st_size > 0
    return False
