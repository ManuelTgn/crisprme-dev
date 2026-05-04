"""
dna_alphabet.py
--------
Nucleotide alphabet definitions, IUPAC encoding tables, and sequence
utility functions shared across the CRISPRme2 package.

Constants
~~~~~~~~~
The module exposes five public constants that downstream code should import
rather than re-define:

- :data:`DNA`          — ordered DNA alphabet including ``N``
- :data:`RNA`          — ordered RNA alphabet
- :data:`IUPAC`        — complete IUPAC ambiguity alphabet (DNA + degenerate)
- :data:`RC`           — reverse-complement mapping (upper and lower case)
- :data:`IUPACTABLE`   — IUPAC character → constituent bases
- :data:`IUPAC_ENCODER` — sorted constituent bases → IUPAC character

Functions
~~~~~~~~~
- :func:`reverse_complement` — reverse-complement a nucleotide string.
- :func:`dna2rna`            — replace ``T``/``t`` with ``U``/``u`` in a string.
- :func:`dna2rna_nt`         — convert a single base character DNA → RNA.
"""

from itertools import permutations
from typing import Dict, List

import os

from .crisprme2_error import Crisprme2ReverseComplementError, Crisprme2DnaRnaError
from .logger import CrisprmeLoggers


# ------------------------------------------------------------------------------
# alphabet constants
# ------------------------------------------------------------------------------

# Ordered DNA alphabet.  ``N`` is included as the fully-ambiguous base.
DNA: List[str] = ["A", "C", "G", "T", "N"]

# Ordered RNA alphabet.  Uracil replaces thymine relative to :data:`DNA`.
RNA: List[str] = ["A", "C", "G", "U"]

# Complete IUPAC nucleotide alphabet: unambiguous DNA bases plus all
# standard degenerate codes (R, Y, S, W, K, M, B, D, H, V).
# ``N`` is the fully-ambiguous wildcard (matches A, C, G, T).
IUPAC: List[str] = DNA + ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]


# ---------------------------------------------------------------------------
# reverse-complement table
# ---------------------------------------------------------------------------

# Reverse-complement mapping for both upper- and lower-case IUPAC characters.
#
# Watson–Crick pairs (A↔T, C↔G) and IUPAC degenerate complements are
# included.  ``U`` is treated as equivalent to ``T`` on the complement
# strand (``U`` → ``A``).  Case is preserved: an upper-case input maps to
# an upper-case complement, and vice-versa.
#
# Used by :func:`reverse_complement`.
RC: Dict[str, str] = {
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
    # lower-case equivalents
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


# ---------------------------------------------------------------------------
# IUPAC encoding tables
# ---------------------------------------------------------------------------

# Mapping from IUPAC ambiguity character to the set of constituent bases
# it represents.
#
# For example, ``IUPACTABLE["R"] == "AG"`` because ``R`` represents either
# purine (A or G).  Single-base entries (A, C, G, T) map to themselves.
#
# This table is the authoritative definition used to derive
# :data:`IUPAC_ENCODER`.
IUPACTABLE: Dict[str, str] = {
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

# Reverse mapping of :data:`IUPACTABLE`: any permutation of a constituent
# base string maps back to its IUPAC character.
#
# Built at import time from :data:`IUPACTABLE` by enumerating all
# permutations of each value string.  For example, both ``"AG"`` and
# ``"GA"`` map to ``"R"``.
#
# Typical use: encode a sorted or unsorted combination of observed bases
# (e.g. from a VCF record) into a single IUPAC character::
#
#     bases = "".join(sorted({"A", "G"}))  # "AG"
#     iupac_char = IUPAC_ENCODER[bases]    # "R"
IUPAC_ENCODER: Dict[str, str] = {
    perm: k
    for k, v in IUPACTABLE.items()
    for perm in {"".join(p) for p in permutations(v)}
}


# ------------------------------------------------------------------------------
# sequence utility functions
# ------------------------------------------------------------------------------


def reverse_complement(sequence: str, loggers: CrisprmeLoggers) -> str:
    """
    Return the reverse complement of a nucleotide string.

    Processes the input right-to-left and maps each character through
    :data:`RC`.  Both upper- and lower-case IUPAC characters are supported;
    case of each base is preserved in the output (upper maps to upper,
    lower to lower).

    Parameters
    ----------
    sequence : str
        Input nucleotide string.  May contain any character present in
        :data:`RC` (standard DNA/RNA bases and IUPAC degenerate codes,
        upper or lower case).
    loggers : CrisprmeLoggers
        Shared logger bundle used for error reporting.

    Returns
    -------
    str
        The reverse complement of *sequence*.

    Raises
    ------
    Crisprme2ReverseComplementError
        If *sequence* contains a character not present in :data:`RC`.

    Examples
    --------
    ::

        rc = reverse_complement("ACGT", loggers)  # "ACGT"
        rc = reverse_complement("GAATTC", loggers) # "GAATTC"  (EcoRI palindrome)
        rc = reverse_complement("NGG", loggers)    # "CCN"
    """
    try:
        return "".join([RC[nt] for nt in sequence[::-1]])
    except (KeyError, Exception) as e:
        loggers.errorlog.log_raise_exception(
            f"Failed reverse complement on {sequence}: {e}",
            os.EX_DATAERR,
            Crisprme2ReverseComplementError,
        )


def dna2rna_nt(base: str) -> str:
    """
    Convert a single DNA base character to its RNA equivalent.

    Replaces ``T`` (thymine) with ``U`` (uracil); all other bases are
    returned uppercased and unchanged.  This is a pure Python function with
    no logging dependency, intended for use in tight inner loops or table
    construction (e.g. :func:`~crisprme2.scores.cfd_score.cfdscore._build_mm_lookup`).

    Parameters
    ----------
    base : str
        A single character representing a DNA nucleotide (case-insensitive).

    Returns
    -------
    str
        Upper-case RNA equivalent of *base* (``"T"`` -> ``"U"``).

    Examples
    --------
    ::

        dna2rna_nt("T")  # "U"
        dna2rna_nt("t")  # "U"
        dna2rna_nt("A")  # "A"
        dna2rna_nt("G")  # "G"
    """
    return "U" if base.upper() in {"T"} else base.upper()


def dna2rna(sequence: str, loggers: CrisprmeLoggers) -> str:
    """
    Convert a DNA sequence string to RNA by replacing thymine with uracil.

    Replaces all occurrences of ``T`` with ``U`` and ``t`` with ``u``,
    preserving case for all other characters.  The conversion is applied
    as a direct string substitution without iterating over individual
    characters, making it efficient for long sequences.

    Parameters
    ----------
    sequence : str
        Input DNA sequence string (may be mixed case).
    loggers : CrisprmeLoggers
        Shared logger bundle used for error reporting.

    Returns
    -------
    str
        RNA equivalent of *sequence* with ``T``/``t`` replaced by ``U``/``u``.

    Raises
    ------
    Crisprme2DnaRnaError
        If the string substitution fails unexpectedly.

    Examples
    --------
    ::

        dna2rna("ACGT", loggers)   # "ACGU"
        dna2rna("acgt", loggers)   # "acgu"
        dna2rna("ACGN", loggers)   # "ACGN"
    """
    try:
        return sequence.replace("T", "U").replace("t", "u")
    except ValueError as e:
        loggers.errorlog.log_raise_exception(
            f"Failed translating DNA to RNA on {sequence}: {e}",
            os.EX_DATAERR,
            Crisprme2DnaRnaError,
        )
