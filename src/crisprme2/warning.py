""" """

from colorama import Fore

import sys


def warning(message: str) -> None:
    """Display a warning message.

    Prints a formatted warning message to standard error.

    Args:
        message: The warning message to display.
    """
    sys.stderr.write(f"{Fore.YELLOW}WARNING: {message}.{Fore.RESET}\n")
