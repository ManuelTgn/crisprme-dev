""" """

from typing import Tuple

import os


# fasta file extensions
FASTAEXTENSIONS = valid_extensions = {"fasta", "fa", "fna", "ffn", "faa", "frn", "fas"}

# bgzip/gzip suffixes recognised on top of a FASTA extension (e.g. ``.fa.gz``)
FASTA_COMPRESSED_SUFFIXES = {"gz", "bgz"}

# fai index file extension
FAI = "fai"

# gzi (bgzip) index file extension
GZI = "gzi"


def find_fai_index(fname: str) -> bool:
    # avoid unexpected crashes due to file location
    fai_index = f"{os.path.abspath(fname)}.{FAI}"
    if os.path.exists(fai_index):  # index must be a non empty file
        return os.path.isfile(fai_index) and os.stat(fai_index).st_size > 0
    return False


def find_gzi_index(fname: str) -> bool:
    # a bgzipped FASTA needs a .gzi companion alongside its .fai
    gzi_index = f"{os.path.abspath(fname)}.{GZI}"
    if os.path.exists(gzi_index):  # index must be a non empty file
        return os.path.isfile(gzi_index) and os.stat(gzi_index).st_size > 0
    return False


def fasta_extension(filepath: str) -> Tuple[str, bool]:
    name = os.path.basename(filepath).lower()
    root, ext = os.path.splitext(name)
    compressed = ext.lstrip(".") in FASTA_COMPRESSED_SUFFIXES
    if compressed:  # peel the compression suffix, inspect the real fasta ext
        _, ext = os.path.splitext(root)
    return ext.lstrip("."), compressed
