"""
This module provides functions for encoding nucleotide sequences into bitset
representations using IUPAC codes.

It enables efficient sequence matching by converting nucleotides and ambiguous
codes into Bitset objects for downstream analysis.
"""

from .utils import IUPAC, print_verbosity
from .crisprme2_error import Crisprme2IupacTableError
from .logger import CrisprmeLoggers
from .bitset import Bitset, SIZE

from typing import List
from time import time

import os


def _encoder(nt: str, position: int, loggers: CrisprmeLoggers) -> Bitset:
    bitset = Bitset(SIZE, loggers)  # 4 - bits encoder
    if nt == IUPAC[0]:  # A - 0001
        bitset.set(0)
    elif nt == IUPAC[1]:  # C - 0010
        bitset.set(1)
    elif nt == IUPAC[2]:  # G - 0100
        bitset.set(2)
    elif nt == IUPAC[3]:  # T - 1000
        bitset.set(3)
    elif nt == IUPAC[4]:  # N - 1111 --> any
        bitset.set_bits("1111")
    elif nt == IUPAC[5]:  # R - 0101 G or A
        bitset.set_bits("0101")
    elif nt == IUPAC[6]:  # Y - 1010 C or T
        bitset.set_bits("1010")
    elif nt == IUPAC[7]:  # S - 0110 C or G
        bitset.set_bits("0110")
    elif nt == IUPAC[8]:  # W  - 1001 A or T
        bitset.set_bits("1001")
    elif nt == IUPAC[9]:  # K - 1100 G or T
        bitset.set_bits("1100")
    elif nt == IUPAC[10]:  # M - 0011 A or C
        bitset.set_bits("0011")
    elif nt == IUPAC[11]:  # B - 1110 --> not A (T or G or C)
        bitset.set_bits("1110")
    elif nt == IUPAC[12]:  # D - 1101 --> not C (A or G or T)
        bitset.set_bits("1101")
    elif nt == IUPAC[13]:  # H - 1011 --> not G (A or C or T)
        bitset.set_bits("1011")
    elif nt == IUPAC[14]:  # V - 0111 --> not T (A or C or G)
        bitset.set_bits("0111")
    else:  # default case
        loggers.errorlog.log_raise_exception(
            f"The nucleotide {nt} at {position} is not a IUPAC character",
            os.EX_DATAERR,
            Crisprme2IupacTableError,
        )
    return bitset


def encode(sequence: str, loggers: CrisprmeLoggers) -> List[Bitset]:
    # encode sequence in bits for efficient matching
    bits = [_encoder(nt.upper(), i, loggers) for i, nt in enumerate(sequence)]
    assert len(bits) == len(sequence)
    return bits
