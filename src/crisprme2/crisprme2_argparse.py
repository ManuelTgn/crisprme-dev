""" """

from .utils import DNA, IUPAC, VCFEXTENSIONS
from .version import __version__
from .fasta import FASTAEXTENSIONS
from .crisprme2 import COMMAND

from argparse import (
    SUPPRESS,
    ArgumentParser,
    HelpFormatter,
    Action,
    _MutuallyExclusiveGroup,
    Namespace,
)
from typing import Iterable, Optional, TypeVar, Tuple, Dict, NoReturn, List, Set
from colorama import Fore
from glob import glob

import multiprocessing
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
        self._initialize_args()  # initialize input args

    def _check_consistency(self):  # sourcery skip: low-code-quality
        """Check the consistency and validity of parsed input arguments.

        Validates the existence, type, and content of input files and directories,
        and sets the list of VCF files found in the specified directory.

        Returns:
            None
        """
        # input fasta files folder
        _validate_directory(
            self._args.genome_dir,
            self._parser,
            f"Cannot find input FASTA folder {self._args.genome_dir}",
        )
        # input vcf files folder
        if self._args.vcf:  # if no input vcf, skip
            _validate_directory(
                self._args.vcf, self._parser, f"Cannot find VCF folder {self._args.vcf}"
            )
        # input guide
        if self._args.fasta_guide:
            _validate_file(
                self._args.fasta_guide,
                self._parser,
                f"Cannot find input guide FASTA {self._args.fasta_guide}",
            )
        elif self._args.bed_guide:
            _validate_file(
                self._args.bed_guide,
                self._parser,
                f"Cannot find input guide BED {self._args.bed_guide}",
            )
        # output folder
        parent_folder = os.path.dirname(self._args.outdir)
        _validate_directory(
            parent_folder, self._parser, f"Cannot find parent folder {parent_folder}"
        )
        # threads number
        _validate_threads(self._args.threads, self._parser)

    def _initialize_args(self) -> None:
        # retreive fasta files in input folder
        self._fastas = _retrieve_files(
            self._args.genome_dir,
            FASTAEXTENSIONS,
            self._parser,
            f"No FASTA file found in {self._args.genome_dir}",
        )
        if self._args.vcf:  # retreive vcf files in input folder
            self._vcfs = _retrieve_files(
                self._args.vcf,
                VCFEXTENSIONS,
                self._parser,
                f"No VCF file found in {self._args.vcf}",
            )
        # retrieve input guide sequence or files
        self._guide, self._guidefasta, self._guidebed = _initialize_guides(
            self._args.guide, self._args.fasta_guide, self._args.bed_guide, self._parser
        )
        # retrieve pam sequence
        self._pam = _initialize_pam(self._args.pam, self._parser)
        # retrieve output folder
        self._outdir = _initialize_outputdir(self._args.outdir)
        # retrieve number of threads
        self._threads = _initialize_threads(self._args.threads)
        # initialize number of mismatches
        self._mm = _initialize_mm(self._args.mm, self._parser)
        # initialize number of DNA/RNA bulges
        self._bdna = _initialize_bdna(self._args.bdna, self._parser)
        self._brna = _initialize_brna(self._args.brna, self._parser)

    @property
    def fastas(self) -> List[str]:
        return self._fastas

    @property
    def vcfs(self) -> List[str]:
        if hasattr(self, "_vcfs"):
            return self._vcfs
        return []

    @property
    def guide(self) -> Optional[str]:
        return self._guide

    @property
    def fasta_guide(self) -> Optional[str]:
        return self._guidefasta

    @property
    def bed_guide(self) -> Optional[str]:
        return self._guidebed

    @property
    def pam(self) -> str:
        return self._args.pam

    @property
    def mm(self) -> int:
        return self._mm

    @property
    def bdna(self) -> int:
        return self._bdna

    @property
    def brna(self) -> int:
        return self._brna

    @property
    def right(self) -> bool:
        return self._args.right

    @property
    def outdir(self) -> str:
        return self._outdir

    @property
    def threads(self) -> int:
        return self._threads


def _validate_directory(
    dirname: str, parser: Crisprme2ArgumentParser, errmsg: str
) -> None:
    if not os.path.exists(dirname) or not os.path.isdir(dirname):  # folder exists?
        parser.error(errmsg)  # print error message to stderr


def _validate_file(fname: str, parser: Crisprme2ArgumentParser, errmsg: str) -> None:
    if not os.path.exists(fname) or not os.path.isfile(fname):  # file exists?
        parser.error(errmsg)  # print error message to stderr


def _validate_threads(threads, parser: Crisprme2ArgumentParser) -> None:
    max_threads = multiprocessing.cpu_count()
    if threads < 0 or threads > max_threads:
        parser.error(
            f"Forbidden number of threads provided ({threads}). Max number of "
            f"available cores: {max_threads}"
        )


def _retrieve_files(
    dirname: str, extensions: Set[str], parser: Crisprme2ArgumentParser, errmsg: str
) -> List[str]:
    fnames = []  # retrieved files list
    for ext in extensions:  # check for each input extension
        fnames.extend(glob(os.path.join(dirname, f"*.{ext}")))
    if not fnames:  # no file found with extensions in folder
        parser.error(errmsg)  # throw error
    fnames = [os.path.abspath(f) for f in fnames]  # avoid ambigous file locations
    return fnames


def _initialize_guides(
    guideseq: str, fasta_guide: str, bed_guide: str, parser: Crisprme2ArgumentParser
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    if guideseq and any(nt.upper() not in DNA[:-1] for nt in guideseq):
        parser.error(f"Invalid guide sequence: {guideseq}")
    guide = guideseq if guideseq else None
    guidefasta = fasta_guide if fasta_guide else None
    guidebed = bed_guide if bed_guide else None
    return guide, guidefasta, guidebed


def _initialize_pam(pam: str, parser: Crisprme2ArgumentParser) -> str:
    if any(nt.upper() not in IUPAC for nt in pam):
        parser.error(f"Invalid PAM sequence {pam}")
    return pam


def _initialize_outputdir(outdir: str) -> str:
    if not os.path.exists(outdir) or not os.path.isdir(outdir):
        os.makedirs(outdir)
    return os.path.abspath(outdir)


def _initialize_threads(threads: int) -> int:
    max_threads = multiprocessing.cpu_count()
    return max_threads if threads == 0 else threads


def _initialize_mm(mm: int, parser: Crisprme2ArgumentParser) -> int:
    if mm < 0:
        parser.error(f"Invalid number of mismatches selected ({mm})")
    return mm


def _initialize_bdna(bdna: int, parser: Crisprme2ArgumentParser) -> int:
    if bdna < 0:
        parser.error(f"Invalid number of DNA bulges selected ({bdna})")
    return bdna


def _initialize_brna(brna: int, parser: Crisprme2ArgumentParser) -> int:
    if brna < 0:
        parser.error(f"Invalid number of RNA bulges selected ({brna})")
    return brna
