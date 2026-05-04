""" """

from itertools import permutations

import os

from .crisprme2_error import Crisprme2ReverseComplementError, Crisprme2DnaRnaError
from .logger import CrisprmeLoggers

# dna alphabet
DNA = ["A", "C", "G", "T", "N"]

# rna alphabet
RNA = ["A", "C", "G", "U"]

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


def reverse_complement(sequence: str, loggers: CrisprmeLoggers) -> str:
    try:
        return "".join([RC[nt] for nt in sequence[::-1]])
    except (KeyError, Exception) as e:
        loggers.errorlog.log_raise_exception(
            f"Failed reverse complement on {sequence}: {e}",
            os.EX_DATAERR,
            Crisprme2ReverseComplementError,
        )


def dna2rna_nt(base: str) -> str:
    """Convert a single DNA base character to its RNA equivalent (T → U)"""
    return "U" if base.upper() in {"T"} else base.upper()


def dna2rna(sequence: str, loggers: CrisprmeLoggers) -> str:
    try:
        return sequence.replace("T", "U").replace("t", "u")
    except ValueError as e:
        loggers.errorlog.log_raise_exception(
            f"Failed translating DNA to RNA on {sequence}: {e}",
            os.EX_DATAERR,
            Crisprme2DnaRnaError,
        )
