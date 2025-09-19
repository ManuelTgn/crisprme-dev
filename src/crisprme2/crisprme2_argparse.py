"""Argument parsing and validation for the CRISPR-HAWK command-line interface.

This module defines custom argument parsers and input argument handler classes for
the CRISPR-HAWK tool, supporting search, VCF conversion, data preparation, and
CRISPRitz configuration workflows. It ensures input consistency, provides helpful
error messages, and exposes validated arguments as convenient properties.
"""

from .utils import COMMAND
from .crisprme2_version import __version__

from argparse import (
    SUPPRESS,
    ArgumentParser,
    HelpFormatter,
    Action,
    _MutuallyExclusiveGroup,
    Namespace,
)
from typing import Iterable, Optional, TypeVar, Tuple, Dict, NoReturn
from colorama import Fore
from glob import glob

import sys
import os

# define abstract generic types for typing
_D = TypeVar("_D")
_V = TypeVar("_V")


class Crisprme2ArgumentParser(ArgumentParser):
    """Custom argument parser for CRISPRme2 command-line interface.

    This class extends argparse.ArgumentParser to provide custom help formatting,
    error handling, and version display for the CRISPRme2 tool.

    Attributes:
        usage (str): The usage string for the parser, with version information.
        formatter_class (type): The custom help formatter class.
    """

    class Crisprme2HelpFormatter(HelpFormatter):
        """Custom help formatter for CRISPRme2 argument parser.

        This formatter customizes the usage message display for the help output.

        Attributes:
            None
        """

        def add_usage(  # type: ignore
            self,
            usage: str,
            actions: Iterable[Action],
            groups: Iterable[_MutuallyExclusiveGroup],
            prefix: Optional[str] = None,
        ) -> None:
            """Add a usage message to the help output.

            Displays the usage description unless suppressed.

            Args:
                usage (str): The usage string to display.
                actions (Iterable[Action]): The actions associated with the parser.
                groups (Iterable[_MutuallyExclusiveGroup]): Mutually exclusive
                    groups.
                prefix (Optional[str]): Optional prefix for the usage message.
            """
            # add usage description for help only if the set action is not to
            # suppress the display of the help formatter
            if usage != SUPPRESS:
                args = (usage, actions, groups, "")
                self._add_item(self._format_usage, args)  # initialize the formatter

    def __init__(self, *args: Tuple[_D], **kwargs: Dict[_D, _V]) -> None:
        """Initialize the CRISPRme2 argument parser.

        Sets up the parser with a custom help formatter and version display.

        Args:
            *args: Positional arguments for ArgumentParser.
            **kwargs: Keyword arguments for ArgumentParser.
        """
        # set custom help formatter defined as
        kwargs["formatter_class"] = self.Crisprme2HelpFormatter  # type: ignore
        # replace the default version display in usage help with a custom
        # version display formatter
        if "usage" in kwargs:
            kwargs["usage"] = kwargs["usage"].replace("{version}", __version__)  # type: ignore
        # initialize argument parser object with input parameters for
        # usage display
        super().__init__(*args, **kwargs)  # type: ignore

    def error(self, error: str) -> NoReturn:  # type: ignore
        """Display an error message and exit.

        Shows the error in red and suggests running the help command.

        Args:
            error (str): The error message to display.

        Raises:
            SystemExit: Exits the program with a usage error code.
        """
        # display error messages raised by argparse in red
        errormsg = (
            f"{Fore.RED}\nERROR: {error}.{Fore.RESET}"
            + f"\n\nRun {COMMAND} -h for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error

    def error_noargs(self) -> None:
        """Display help and exit when no arguments are provided.

        Prints the help message and exits with a no input code.

        Raises:
            SystemExit: Exits the program with a no input error code.
        """
        self.print_help()  # if no input argument, print help
        sys.exit(os.EX_NOINPUT)  # exit with no input code


class Crisprme2SearchInputArgs:
    """Handles and validates parsed command-line arguments for CRISPRme2.

    This class checks the consistency of input arguments and provides convenient
    access to validated argument values as properties.

    Attributes:
        _args (Namespace): The parsed arguments namespace.
        _parser (Crisprme2ArgumentParser): The argument parser instance.
    """

    def __init__(self, args: Namespace, parser: Crisprme2ArgumentParser) -> None:
        """Initialize Crisprme2SearchInputArgs with parsed arguments and parser.

        Stores the parsed arguments and parser, then checks argument consistency.

        Args:
            args (Namespace): The parsed arguments namespace.
            parser (Crisprme2ArgumentParser): The argument parser instance.
        """
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _check_consistency(self):  # sourcery skip: low-code-quality
        """Check the consistency and validity of parsed input arguments.

        Validates the existence, type, and content of input files and directories,
        and sets the list of VCF files found in the specified directory.

        Returns:
            None
        """
        # fasta file
        if not os.path.exists(self._args.genome_dir) or not os.path.isdir(self._args.genome_dir):
            self._parser.error(f"Cannot find input FASTA folder {self._args.fasta}")
        self._fastas = glob(os.path.join(self._args.genome_dir, "*.fa")) + glob(
            os.path.join(self._args.genome_dir, "*.fasta")
        )
        if not self._fastas:
            self._parser.error(f"No FASTA file found in {self._args.genome_dir}")
        # vcf folder
        if self._args.vcf and (not os.path.isdir(self._args.vcf)):
            self._parser.error(f"Cannot find VCF folder {self._args.vcf}")
        self._vcfs = glob(os.path.join(self._args.vcf, "*.vcf.gz"))
        if self._args.vcf and not self._vcfs:
            self._parser.error(f"No VCF file found in {self._args.vcf}")
        # output folder
        parent_folder = os.path.dirname(self._args.outdir)
        if not os.path.exists(self._args.outdir) or not os.path.isdir(
            self._args.outdir
        ):
            if not os.path.exists(parent_folder) or not os.path.isdir(parent_folder):
                self._parser.error(f"Cannot find output folder {self._args.outdir}")
            os.makedirs(self._args.outdir)
    