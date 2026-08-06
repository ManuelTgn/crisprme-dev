""" """

from pathlib import Path
from pysam import faidx, FastaFile
from pysam.utils import SamtoolsError
from typing import Dict, List, Optional

import contextlib
import os


from .crisprme2_error import Crisprme2FastaError, Crisprme2SequenceError
from .fasta_utils import (
    fasta_extension,
    find_fai_index,
    find_gzi_index,
    FAI,
    FASTAEXTENSIONS,
)
from .logger import CrisprmeLoggers
from .sequence import ContigSequence
from .warning import warning


class Fasta:

    def __init__(self, filepath: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._filepath = filepath  # fasta filename
        self._validate_file()  # validate fasta file structure
        self._index = self._search_index()  # fai index
        self._fasta_handle: Optional[FastaFile] = None
        self._is_open = False
        self._init_contig_length()  # initialize contig name(s) and length(s)

    def _validate_file(self) -> None:
        # check file extension; bgzipped FASTA (.fa.gz) is accepted. The
        # compression flag is stored for .gzi handling in _search_index
        ext, self._compressed = fasta_extension(self._filepath)
        if ext not in FASTAEXTENSIONS:
            self._loggers.errorlog.log_raise_exception(
                f"File {self._filepath} does not have a standard FASTA extension",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )

    def _index_fasta(self, pytest: bool = False) -> str:
        if hasattr(self, "_index") and (self._index and not pytest):
            warning("FASTA index already present, forcing update")
        try:  # create index in the same folder as the input fasta
            self._loggers.verboselog.debug(
                f"Creating index for FASTA: {self._filepath}"
            )
            faidx(str(self._filepath))  # bgzipped input also yields a .gzi
        except (OSError, Exception) as e:
            self._loggers.errorlog.log_exception(
                f"Failed indexing FASTA {self._filepath}: {e} "
                "(compressed FASTA must be bgzip-compressed, not plain gzip)",
                os.EX_DATAERR,
            )
        assert find_fai_index(str(self._filepath))  # now should be available
        if self._compressed:  # bgzipped -> .gzi must sit alongside the .fai
            assert find_gzi_index(str(self._filepath))
        return f"{self._filepath}.{FAI}"

    def _search_index(self) -> Path:
        # a plain FASTA needs only a .fai; a bgzipped FASTA needs a .gzi too.
        # (Re)compute the index when either companion is missing
        fai_present = find_fai_index(str(self._filepath))
        gzi_present = (not self._compressed) or find_gzi_index(str(self._filepath))
        if fai_present and gzi_present:  # indexes present, store the .fai path
            return Path(f"{self._filepath}.{FAI}")
        # index not found -> compute it de novo and store it in the same folder
        # as the input fasta
        self._loggers.verboselog.debug(f"FASTA index not found for {self._filepath}")
        return Path(self._index_fasta())

    def _init_contig_length(self) -> None:
        self.open()  # manually open fasta file
        assert self._fasta_handle  # should be open
        self._ncontigs = len(self._fasta_handle.references)
        # store name and length for every contig in the fasta, so both
        # single- and multi-contig fasta files are supported uniformly
        self._contigs = list(self._fasta_handle.references)
        self._lengths: Dict[str, int] = dict(
            zip(self._fasta_handle.references, self._fasta_handle.lengths)
        )
        self.close()  # manually close fasta file

    def open(self) -> "Fasta":
        if self._is_open:
            self._loggers.errorlog.log_raise_exception(
                f"FASTA file {self._filepath} is already open",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        try:  # open fasta, assumes that index is already available
            self._fasta_handle = FastaFile(str(self._filepath))
            self._is_open = True
        except (OSError, Exception) as e:
            self._loggers.errorlog.log_exception(
                f"Failed to open FASTA file {self._filepath}: {str(e)}", os.EX_IOERR
            )
        return self

    def close(self) -> None:
        if self._fasta_handle is not None:
            self._fasta_handle.close()
            self._is_open = False

    def __enter__(self) -> "Fasta":
        return self.open()

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def read(self) -> bytearray:
        """Return the raw sequence bytes in the fasta file, with every
        header ('>') line stripped out and all contigs' sequences
        concatenated together. For per-contig access, use `read_contigs`.
        """
        sequence = bytearray()
        with open(self._filepath, mode="rb") as fin:
            for line in fin:
                if line.startswith(b">"):
                    continue  # skip header line(s)
                sequence.extend(line.rstrip(b"\r\n"))
        return sequence

    def read_contigs(self) -> Dict[str, bytearray]:
        sequences: Dict[str, bytearray] = {}
        current: Optional[str] = None
        with open(self._filepath, mode="rb") as fin:
            for line in fin:
                if line.startswith(b">"):
                    current = line[1:].strip().split()[0].decode()
                    sequences[current] = bytearray()
                elif current is not None:
                    sequences[current].extend(line.rstrip(b"\r\n"))
        return sequences

    def fetch(
        self, reference: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> ContigSequence:
        if not self._is_open or self._fasta_handle is None:
            self._loggers.errorlog.log_raise_exception(
                "FASTA file must be opened before fetching",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        assert self._fasta_handle  # must not be none
        try:
            if start is None and end is None:  # access string by contig name
                return ContigSequence(
                    self._fasta_handle.fetch(reference),
                    reference,
                    0,
                    self._lengths[reference],
                    self._loggers,
                )
            elif start is not None and end is not None:
                if start < 0 or end < start:
                    self._loggers.errorlog.log_raise_exception(
                        f"Invalid coordinates: start={start}, end={end}",
                        os.EX_DATAERR,
                        Crisprme2SequenceError,
                    )
                return ContigSequence(
                    self._fasta_handle.fetch(reference, start, end),
                    reference,
                    start,
                    end,
                    self._loggers,
                )
            else:
                self._loggers.errorlog.log_raise_exception(
                    "Both start and end must be specified or both None",
                    os.EX_DATAERR,
                    Crisprme2SequenceError,
                )
        except KeyError:
            self._loggers.errorlog.log_raise_exception(
                f"Reference '{reference}' not found in FASTA file",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        except Exception as e:
            self._loggers.errorlog.log_exception(
                f"Error fetching sequence: {str(e)}", os.EX_DATAERR
            )

    def __contains__(self, reference: str) -> bool:
        return (
            reference in self._fasta_handle.references
            if self._is_open and self._fasta_handle
            else False
        )

    def __repr__(self) -> str:
        status = "open" if self._is_open else "closed"
        sequences = (
            self._contigs[0] if self._ncontigs == 1 else f"{self._ncontigs} contigs"
        )
        return f"<{self.__class__.__name__} object; sequences={sequences}, status={status}>"

    def __del__(self):
        if self._is_open:
            self.close()

    @property
    def contig(self) -> str:
        """Name of the contig, 'chr'-prefixed.

        Only valid for single-contig fasta files; use `contigs` for
        fasta files with multiple contigs.
        """
        if self._ncontigs != 1:
            self._loggers.errorlog.log_raise_exception(
                f"FASTA file {self._filepath} contains {self._ncontigs} contigs; "
                "use 'contigs' instead of 'contig'",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        return self._contigs[0]

    @property
    def contigs(self) -> List[str]:
        # return every contig name in the fasta file, 'chr'-prefixed
        return self._contigs

    @property
    def length(self) -> int:
        """Length of the contig.

        Only valid for single-contig fasta files; use `lengths` or
        `get_length(reference)` for fasta files with multiple contigs.
        """
        if self._ncontigs != 1:
            self._loggers.errorlog.log_raise_exception(
                f"FASTA file {self._filepath} contains {self._ncontigs} contigs; "
                "use 'lengths' or 'get_length(reference)' instead of 'length'",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        return self._lengths[self._contigs[0]]

    @property
    def lengths(self) -> Dict[str, int]:
        # return a contig name -> length mapping, for every contig
        return dict(self._lengths)

    def get_length(self, reference: str) -> int:
        # return the length of one specific contig/reference
        if reference not in self._lengths:
            self._loggers.errorlog.log_raise_exception(
                f"Reference '{reference}' not found in FASTA file",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        return self._lengths[reference]

    @property
    def nreferences(self) -> int:
        return self._ncontigs  # return the number of contigs in input fasta


class GuideFasta(Fasta):

    def __init__(self, filepath: str, loggers: CrisprmeLoggers) -> None:
        super().__init__(filepath, loggers)

    def _read_guides(self) -> None:
        f = FastaFile(str(self._filepath), filepath_index=str(self._index))
        gnames = f.references  # retrieve guides seqnames (fasta headers)
        try:  # extract guide sequence
            self._guides = list({f.fetch(gname) for gname in gnames})
        except (SamtoolsError, Exception):
            self._loggers.errorlog.log_exception(
                f"Failed parsing guides from {self._filepath}", os.EX_DATAERR
            )

    @property
    def guides(self) -> List[str]:
        return self._guides


def read_fasta_files(
    fasta_files: List[str], loggers: CrisprmeLoggers
) -> Dict[str, Fasta]:
    fastas: Dict[str, Fasta] = {}  # fasta-contig map
    for fasta_file in fasta_files:
        loggers.verboselog.debug(f"Create FASTA object {fasta_file}")
        try:  # validates + ensures index + contig/length
            fasta = Fasta(fasta_file, loggers)
            contigs = fasta.contigs
        except Exception:  # Fasta() might have opened internally -> close
            with contextlib.suppress(Exception):
                fasta.close()  # type: ignore[name-defined]
            loggers.errorlog.log_raise_exception(
                f"Failed FASTA object creation: {fasta_file}", os.EX_IOERR, IOError
            )
        for contig in contigs:
            if contig in fastas:
                loggers.errorlog.log_raise_exception(
                    f"Multiple FASTA files with contig {contig}",
                    os.EX_DATAERR,
                    Crisprme2FastaError,
                )
            fastas[contig] = fasta
        loggers.verboselog.debug(f"Successfully FASTA object created: {fasta_file}")
    return fastas
