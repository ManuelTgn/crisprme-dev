""" """

from .logger import CrisprmeLoggers

ENCODING = {'A': 0x0, 'C': 0x1, 'G': 0x2, 'T': 0x3, 'R': 0x4, 'Y': 0x5, 'S': 0x6, 'W': 0x7, 'K': 0x8, 'M': 0x9, 'B': 0xA, 'D': 0xB, 'H': 0xC, 'V': 0xD, 'N': 0xE}
DECODING = ["A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "-"]


class BitSequence:
    def __init__(self, sequence: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._length = len(sequence)
        self._data = bytearray(self._length)  # encoder
        self._encode(sequence)  # encode input string (1 byte per nt) 
        

    def _encode(self, sequence: str) -> None:
        for i, nt in enumerate(sequence):
            self._data[i] = ENCODING.get(nt, ENCODING["N"])

    def decode(self) -> str:
        return "".join(DECODING[self._data[i] & 0xF] for i in range(self._length))