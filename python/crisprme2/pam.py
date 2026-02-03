"""Provides the PAM class for representing and encoding Protospacer Adjacent
Motif sequences.

This module defines the PAM class, which validates, and stores PAM sequences and
their reverse complements for efficient sequence matching.
"""

from .crisprme2_error import Crisprme2PamError
from .logger import CrisprmeLoggers
from .encoder import BitSequence
from .utils import reverse_complement

from typing import List
from time import time

import os

# list PAMs for each cas system
CASXPAM = ["TTCN"]
CPF1PAM = [
    "TTN",
    "TTTN",
    "TYCV",
    "TATV",
    "TTTV",
    "TTTR",
    "ATTN",
    "TTTA",
    "TCTA",
    "TCCA",
    "CCCA",
    "YTTV",
    "TTYN",
]
SACAS9PAM = ["NNGRRT", "NNNRRT"]
SPCAS9PAM = ["NGG", "NGA", "NRG", "NGC"]
XCAS9PAM = ["NGK", "NGN", "NNG"]

# list cas systems
CASX = 0
CPF1 = 1
SACAS9 = 2
SPCAS9 = 3
XCAS9 = 4


class PAM:

    def __init__(self, pamseq: str, right: bool, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._sequence = pamseq.upper()  # store pam sequence
        self._reverse_complement()  # store pam reverse complement
        self._assess_cas_system(right)  # assess cas system
        self._encode()  # encode pam in bits (one byte per nt)

    def __len__(self) -> int:
        """Returns the length of the PAM sequence.

        This method allows the PAM object to be used with the built-in len()
        function.

        Returns:
            int: The length of the PAM sequence.
        """
        return len(self._sequence)

    def __eq__(self, value: object) -> bool:
        """Checks equality between this PAM object and another.

        Compares the stored PAM sequence with another PAM object's sequence to
        determine equality.

        Args:
            value: The object to compare with this PAM instance.

        Returns:
            bool: True if the sequences are equal and the object is a PAM instance,
                False otherwise.
        """
        return self._sequence == value.pam if isinstance(value, PAM) else NotImplemented

    def __repr__(self) -> str:
        """Returns a string representation of the PAM object for debugging.

        This method provides a detailed string useful for developers to inspect
        the PAM object.

        Returns:
            str: A string representation of the PAM object.
        """
        return f"<{self.__class__.__name__} object; sequence={self._sequence}>"

    def __str__(self) -> str:
        """Returns the PAM sequence as a string.

        This method allows the PAM object to be converted to its sequence string
        representation.

        Returns:
            str: The PAM sequence.
        """
        return f"{self._sequence}"

    def _assess_cas_system(self, right: bool) -> None:
        self._cas_system = -1  # unknown cas system pam
        if self._sequence in CASXPAM:  # casx system pam
            self._cas_system = CASX
        elif self._sequence in CPF1PAM and right:  # cpf1 cas system pam
            self._cas_system = CPF1
        elif self._sequence in SACAS9PAM:  # sacas9 system pam
            self._cas_system = SACAS9
        elif self._sequence in SPCAS9PAM and not right:  # spcas9 system pam
            self._cas_system = SPCAS9
        elif self._sequence in XCAS9PAM and not right:  # xcas9 pam
            self._cas_system = XCAS9

    def _reverse_complement(self) -> None:
        assert hasattr(self, "_sequence")  # required
        try:  # reverse complement is used to find off-targets on 3'-5'
            self._sequence_rc = reverse_complement(self._sequence)
        except (KeyError, Exception):
            self._loggers.errorlog.log_exception(
                f"Failed reverse complement on PAM {self._sequence}", os.EX_DATAERR
            )

    def _encode(self) -> None:
        # encode forward and reverse pam
        assert hasattr(self, "_sequence") and hasattr(self, "_sequence_rc")
        self._bitsequence = BitSequence(self._sequence, self._loggers)
        self._bitsequence_rc = BitSequence(self._sequence_rc, self._loggers)

    def decode(self, strand: int) -> str:
        if strand not in [0, 1]:  # unknown strand
            self._loggers.errorlog.log_raise_exception(
                "Only 0 (forward) and 1 (reverse) are values allowed for "
                f"strandness, got {strand}",
                os.EX_DATAERR,
                Crisprme2PamError,
            )
        return (
            self._bitsequence.decode() if strand == 0 else self._bitsequence_rc.decode()
        )

    @property
    def pam(self) -> str:
        return self._sequence

    @property
    def pamb(self) -> bytearray:
        return self._bitsequence.data

    @property
    def rc(self) -> str:
        return self._sequence_rc

    @property
    def rcb(self) -> bytearray:
        return self._bitsequence_rc.data

    @property
    def cas_system(self) -> int:
        return self._cas_system


def read_pam(pamseq: str, loggers: CrisprmeLoggers) -> PAM:
    loggers.verboselog.debug(f"Creating PAM object for PAM {pamseq}")
    start = time()
    pam = PAM(pamseq, False, loggers)  # initialize pam object
    loggers.verboselog.debug(
        f"PAM object for PAM {pam} created in {time() - start:.2f}s"
    )
    return pam
