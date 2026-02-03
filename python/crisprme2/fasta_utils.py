""" """

from .crisprme2_error import Crisprme2FastaError
from .fasta import Fasta
from .logger import CrisprmeLoggers

from typing import List, Dict

import os


def read_fasta_files(
    fasta_files: List[str], loggers: CrisprmeLoggers
) -> Dict[str, Fasta]:
    fastas: Dict[str, Fasta] = {}  # fasta-contig map
    for fasta_file in fasta_files:
        loggers.verboselog.debug(f"Create FASTA object {fasta_file}")
        try:
            fasta = Fasta(
                fasta_file, loggers
            )  # validates + ensures index + contig/length
            contig = fasta.contig
        except Exception:  # Fasta() might have opened interbally -> close
            try:
                fasta.close()  # type: ignore[name-defined]
            except Exception:
                pass
            loggers.errorlog.log_raise_exception(
                f"Failed FASTA object creation: {fasta_file}", os.EX_IOERR, IOError
            )
        if contig in fastas:
            loggers.errorlog.log_raise_exception(
                f"Multiple FASTA files with contig {contig}",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        fastas[contig] = fasta
        loggers.verboselog.debug(f"Successfully FASTA object created: {fasta_file}")
    return fastas
