""" """

from .crisprme2_argparse import Crisprme2SearchInputArgs
from .crisprme2_error import Crisprme2GuideError
from .sequence import Sequence
from .logger import CrisprmeLoggers
from .utils import reverse_complement
from .encoder import BitSequence
from .sequence import Sequence
from .fasta import GuideFasta
from .pam import PAM

from typing import Union, List, Optional
from time import time

import sys
import os


class Guide(Sequence):

    def __init__(self, sequence: str, loggers: CrisprmeLoggers):
        super().__init__(sequence, loggers)
        self._reverse_complement()  # compute guide reverse complement

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; sequence={self.sequence}>"

    def __len__(self) -> int:
        return self._length

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, Guide):
            return NotImplemented
        return self._sequence == value.sequence

    def __getitem__(self, idx: Union[int, slice]) -> Union[str, List[str]]:
        if not hasattr(self, "_sequence"):
            self._loggers.errorlog.log_raise_exception(
                f"Missing _sequence attribute on class {self.__class__.__name__}",
                os.EX_DATAERR,
                AttributeError,
            )
        try:
            return self._sequence[idx]
        except IndexError:
            self._loggers.errorlog.log_exception(
                f"Index {idx} out of bounds", os.EX_DATAERR
            )
            sys.exit(os.EX_DATAERR)

    def _reverse_complement(self) -> None:
        assert hasattr(self, "_sequence")  # required
        try:  # reverse complement is used to find off-targets on 3'-5'
            self._sequence_rc = Sequence(
                reverse_complement(self.sequence), self._loggers
            )
        except (KeyError, Exception):
            self._loggers.errorlog.log_exception(
                f"Failed reverse complement on guide {self.sequence}", os.EX_DATAERR
            )

    def _encode(self) -> None:
        # encode forward and reverse pam
        assert hasattr(self, "_sequence") and hasattr(self, "_sequence_rc")
        self._bitsequence = BitSequence(self.sequence, self._loggers)
        self._bitsequence_rc = BitSequence(self.sequence, self._loggers)

    def decode(self, strand: int) -> str:
        if strand not in [0, 1]:  # unknown strand
            self._loggers.errorlog.log_raise_exception(
                "Only 0 (forward) and 1 (reverse) are values allowed for "
                f"strandness, got {strand}", 
                os.EX_DATAERR, 
                Crisprme2GuideError,
            )
        return self._bitsequence.decode() if strand == 0 else self._bitsequence_rc.decode()


    @property
    def pamb(self) -> bytearray:
        return self._bitsequence.data

    @property
    def rc(self) -> Sequence:
        return self._sequence_rc
    
    @property
    def rcb(self) -> bytearray:
        return self._bitsequence_rc.data


class GuidesList:

    def __init__(self, guides: List[Guide], loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        if any(len(guides[0]) != len(g) for g in guides):
            self._loggers.errorlog.log_raise_exception(
                "Found input guides with different length, provide only guides "
                "sharing the same length", 
                os.EX_DATAERR, 
                Crisprme2GuideError,
            )
        self._guides = guides  # guides list

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; num guides={len(self)}>"

    def __str__(self) -> str:
        return "\n".join(guide.sequence for guide in self)

    def __len__(self) -> int:
        return len(self._guides)

    def __iter__(self) -> "GuidesListIterator":
        return GuidesListIterator(self)

    def __getitem__(self, idx: Union[int, slice]) -> Union[Guide, List[Guide]]:
        if not hasattr(self, "_guides"):
            self._loggers.errorlog.log_raise_exception(
                f"Missing _guides attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                AttributeError,
            )
        try:
            return self._guides[idx]
        except IndexError:
            self._loggers.errorlog.log_exception(
                f"Index {idx} out of range", os.EX_DATAERR
            )
            sys.exit(os.EX_DATAERR)

    def extend(self, value: "GuidesList") -> None:
        if not isinstance(value, GuidesList):
            self._loggers.errorlog.log_raise_exception(
                f"Cannot extend {self.__class__.__name__} with objects of type {type(value).__name__}",
                os.EX_DATAERR,
                TypeError,
            )
        self._guides.extend(value.guides)  # extend guides list

    def append(self, guide: Guide) -> None:
        if not isinstance(guide, Guide):
            self._loggers.errorlog.log_raise_exception(
                f"Cannot append to {self.__class__.__name__} objects of type {type(guide).__name__}",
                os.EX_DATAERR,
                TypeError,
            )
        self._guides.append(guide)

    @property
    def guides(self) -> List[Guide]:
        return self._guides


class GuidesListIterator:

    def __init__(self, guides: GuidesList) -> None:
        self._guides = guides  # guides list object to iterate over
        self._index = 0  # iterator index used over the list

    def __next__(self) -> Guide:
        if self._index < len(self._guides):
            result = self._guides[self._index]
            assert isinstance(result, Guide)
            self._index += 1  # go to next position in list
            return result
        raise StopIteration  # stop iteration over regions list


def _read_guide(guide: str, loggers: CrisprmeLoggers) -> GuidesList:
    # single grna as input (--guide)
    loggers.verboselog.debug(f"Reading input guide: {guide}")
    start = time()
    loggers.verboselog.debug(f"Input guide {guide} read in {time() - start:.2f}s")
    return GuidesList([Guide(guide, loggers)], loggers)


def _read_guides_fasta(fasta_guides: str, loggers: CrisprmeLoggers) -> GuidesList:
    # multiple grnas as input in fasta (--sequence)
    loggers.verboselog.debug(f"Reading input guide FASTA: {fasta_guides}")
    start = time()
    gf = GuideFasta(fasta_guides, loggers)  # parse fasta file
    guides = [Guide(guide, loggers) for guide in gf.guides]
    loggers.verboselog.debug(
        f"Read {len(guides)} guides in {fasta_guides} in {time() - start:.2f}s"
    )
    return GuidesList(guides, loggers)


def read_guides(args: Crisprme2SearchInputArgs, loggers: CrisprmeLoggers) -> GuidesList:
    # only one option is allowed
    assert sum(bool(e) for e in [args.guide, args.fasta_guide, args.bed_guide]) == 1
    if args.guide:  # --guide option (single guide)
        return _read_guide(args.guide, loggers)
    elif args.fasta_guide:  # extract guide sequence
        return _read_guides_fasta(args.fasta_guide, loggers)
    elif args.bed_guide:  # --coordinates option (guides extracted via bed)
        pass
    loggers.errorlog.log_raise_exception(
        "Invalid input: no guide input option selected. None of the following "
        "selected: --guide, --sequence, or --coordinates",
        os.EX_DATAERR,
        Crisprme2GuideError,
    )
    sys.exit(os.EX_DATAERR)
