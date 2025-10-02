"""Provides a Bitset class for efficient manipulation of fixed-size binary data.

This module defines the Bitset class, which allows for setting, resetting, testing,
and performing bitwise operations on individual bits.
It is designed for compact storage and manipulation of binary data, with error
handling for invalid operations.
"""

from .logger import CrisprmeLoggers
from .crisprme2_error import Crisprme2BitsetError
from .exception_handlers import exception_handler

import os

SIZE = 4  # bit size


class Bitset(object):
    def __init__(self, size: int, loggers: CrisprmeLoggers) -> None:
        if size < 1:
            loggers.errorlog.log_raise_exception(
                f"Forbidden Bitset size ({size})", os.EX_DATAERR, Crisprme2BitsetError
            )
        self._loggers = loggers  # store loggers
        self._size = size  # number of bits to encode
        self._bits = 0  # bit initialized to 0

    def __str__(self) -> str:
        """Return the bitset as a zero-padded binary string.

        Converts the internal integer representation to a binary string, padded
        to the bitset's size.

        Returns:
            The binary string representation of the bitset.
        """
        # bin(self._bits) converts an integer to a binary string
        # [2:] remove the first two characters from the output of bin(self._bits)
        # .zfill(self._size) add as many zeros as required to reach self._size
        # Example: bits = 5, size =  8 bits  --> 101 --> 00000101 final
        return bin(self._bits)[2:].zfill(self._size)

    def __repr__(self) -> str:
        """Return a string representation of the Bitset object.

        Provides a detailed string showing the class name, binary value, and
        size of the bitset.

        Returns:
            The string representation of the Bitset.
        """
        # value is the binary string, size is the len of the binary string
        return f"<{self.__class__.__name__} object; value={self}, size={self._size}>"

    def __and__(self, bitset: "Bitset") -> "Bitset":
        """Perform a bitwise AND operation with another Bitset.

        Returns a new Bitset that is the result of the bitwise AND between this
        bitset and another of the same size.

        Args:
            bitset: Another Bitset object to perform the AND operation with.

        Returns:
            A new Bitset object representing the result of the AND operation.

        Raises:
            Crisprme2BitsetError: If the two bitsets are not the same size.
        """
        if self._size != bitset.size:
            self._loggers.errorlog.log_raise_exception(
                f"{self.__class__.__name__} objects must have the same size for AND operator",
                os.EX_DATAERR,
                Crisprme2BitsetError,
            )
        result = Bitset(self._size, self._loggers)  # allocate bits for AND result
        result._bits = self._bits & bitset.bits  # perform AND between bitsets
        return result

    @property
    def size(self) -> int:
        return self._size

    @property
    def bits(self) -> int:
        return self._bits

    def set(self, index: int) -> None:
        """Set the bit at the specified index to 1.

        Updates the bit at the given index to 1, raising an error if the index
        is out of bounds.

        Args:
            index: The position of the bit to set.

        Raises:
            Crisprme2BitsetError: If the index is out of bounds.
        """
        if index >= self._size:
            self._loggers.errorlog.log_raise_exception(
                f"Index {index} out of bounds, unable to set bit",
                os.EX_DATAERR,
                Crisprme2BitsetError,
            )
        # bitwise OR operation -> sets 1 at the input position and shifts the
        # 1s to the left (<<)
        # EXAMPLE: if index = 3, 0000 -> 1000
        self._bits |= 1 << index

    def reset(self, index: int) -> None:
        """Reset the bit at the specified index to 0.

        Sets the bit at the given index to 0, raising an error if the index is
        out of bounds.

        Args:
            index: The position of the bit to reset.

        Raises:
            Crisprme2BitsetError: If the index is out of bounds.
        """
        if index >= self._size:
            self._loggers.errorlog.log_raise_exception(
                f"Index {index} out of bounds, unable to reset bit",
                os.EX_DATAERR,
                Crisprme2BitsetError,
            )
        # reset bit at position index
        self._bits &= ~(1 << index)  # ~ is not operator

    def set_bits(self, bits: str) -> None:
        """Set the bits of the bitset using a binary string.

        Updates the bitset to match the provided string of '0's and '1's,
        raising an error for invalid input.

        Args:
            bits: A string representing the bits to set, consisting only of '0'
                and '1'.

        Raises:
            Crisprme2BitsetError: If the input string contains characters other
                than '0' or '1'.
        """
        if any(bit not in "01" for bit in bits):
            self._loggers.errorlog.log_raise_exception(
                f"{bits} is not a bit string",
                os.EX_DATAERR,
                Crisprme2BitsetError,
            )
        bitstring_size = len(bits)
        for i, bit in enumerate(bits):
            if bit == "0":  # force bit reset
                self.reset(bitstring_size - 1 - i)
            else:  # force set bit (bit == 1)
                self.set(bitstring_size - 1 - i)
        assert str(self) == bits

    def to_bool(self) -> bool:
        """Convert the bitset to a boolean value.

        Returns True if any bit in the bitset is set to 1, otherwise returns False.

        Returns:
            True if the bitset contains at least one set bit, False otherwise.
        """
        return bool(self._bits)  # cast bitset to bool

    def test(self, index: int) -> bool:
        """Test whether the bit at the specified index is set to 1.

        Checks if the bit at the given index is set, raising an error if the
        index is out of bounds.

        Args:
            index: The position of the bit to test.

        Returns:
            True if the bit at the specified index is set to 1, False otherwise.

        Raises:
            Crisprme2BitsetError: If the index is out of bounds.
        """
        if index >= self._size:
            self._loggers.errorlog.log_raise_exception(
                f"Index {index} out of bounds, unable to test bit",
                os.EX_DATAERR,
                Crisprme2BitsetError,
            )
        return bool(self._bits & (1 << index))  # test if bit at position index
