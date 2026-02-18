""" """

from .verbosity import VERBOSITYLVL

from colorama import Fore

import sys


def warning(message: str, verbosity: int) -> None:
    """Display a warning message if the verbosity level is sufficient.

    Prints a formatted warning message to standard error if the verbosity
    threshold is met.

    Args:
        message: The warning message to display.
        verbosity: The current verbosity level.
    """
    if verbosity >= VERBOSITYLVL[1]:
        sys.stderr.write(f"{Fore.YELLOW}WARNING: {message}.{Fore.RESET}\n")
    return
