""" """

from colorama import Fore

import sys

# define static variables shared across software modules
TOOLNAME = "CRISPRme2"  # tool name
COMMAND = "crisprme2"  # command line call
# define verbosity levels
VERBOSITYLVL = [0, 1, 2, 3]


def print_verbosity(message: str, verbosity: int, verbosity_threshold: int) -> None:
    """Print a message if the verbosity level meets the threshold.

    Outputs the provided message to standard output if the current verbosity is
    greater than or equal to the specified threshold.

    Args:
        message: The message to print.
        verbosity: The current verbosity level.
        verbosity_threshold: The minimum verbosity level required to print the
            message.
    """
    if verbosity >= verbosity_threshold:
        sys.stdout.write(f"{message}\n")
    return

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
