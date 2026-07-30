""" """

from argparse import Namespace
from glob import glob
from typing import List, Optional

import multiprocessing
import os


from .crisprme2_argparse import Crisprme2ArgumentParser
from .dna_alphabet import DNA, IUPAC


class Crisprme2InputArgs:

    def __init__(self, args: Namespace, parser: Crisprme2ArgumentParser) -> None:
        self._args = args
        self._parser = parser

    def _validate_output_folder(self) -> None:
        outdir = os.path.abspath(self._args.outdir)
        parentdir = os.path.dirname(self._args.outdir) or os.getcwd()
        _check_folder(
            parentdir, self._parser, f"Cannot create output folder {self._args.outdir}"
        )
        os.makedirs(outdir, exist_ok=True)  # create output folder
        self._outdir = outdir

    def _validate_threads(self) -> None:
        """Validate and store the requested thread count.

        Returns
        -------
        None
        """
        self._threads = _validate_threads_num(self._args.threads, self._parser)

    @property
    def outdir(self) -> str:
        """str: Absolute path of the validated output folder."""
        return self._outdir

    @property
    def threads(self) -> int:
        """int: Requested number of threads."""
        return self._args.threads


class Crisprme2SearchInputArgs(Crisprme2InputArgs):

    def __init__(self, args: Namespace, parser: Crisprme2ArgumentParser) -> None:
        """Initialize Crisprme2SearchInputArgs with parsed arguments and parser.

        Stores the parsed arguments and parser, then checks argument consistency.

        Args:
            args (Namespace): The parsed arguments namespace.
            parser (Crisprme2ArgumentParser): The argument parser instance.
        """
        super().__init__(args, parser)
        self._check_consistency()

    def _validate_genome_folder(self) -> None:
        """Validate the genome folder and collect its FASTA files.

        Returns
        -------
        None
        """
        _check_folder(
            self._args.genome,
            self._parser,
            f"Cannot find input genome folder {self._args.genome}",
        )
        self._fastas = glob(os.path.join(self._args.genome, "*.fa")) + glob(
            os.path.join(self._args.genome, "*.fasta")
        )
        _check_retrieved_files(
            self._fastas, self._parser, f"No FASTA file found in {self._args.genome}"
        )

    def _validate_vcf_folder(self) -> None:
        """Validate the VCF folder and collect its ``*.vcf.gz`` files.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            Via the parser's ``error`` method if the folder is missing or
            contains no VCF files.
        """
        _check_folder(
            self._args.vcf, self._parser, f"Cannot find VCF folder {self._args.vcf}"
        )
        self._vcfs = glob(os.path.join(self._args.vcf, "*.vcf.gz"))
        _check_retrieved_files(
            self._vcfs, self._parser, f"No VCF file found in {self._args.vcf}"
        )

    def _validate_guides(self) -> None:
        self._guide, self._guidefasta, self._guidebed = None, None, None
        if self._args.guide:
            if any(nt.upper() not in DNA[:-1] for nt in self._args.guide):
                self._parser.error(f"Invalid guide sequence: {self._args.guide}")
            self._guide = self._args.guide
        if self._args.fasta_guide:
            _check_file(
                self._args.fasta_guide,
                self._parser,
                f"Cannot find input guide FASTA {self._args.fasta_guide}",
            )
            self._guidefasta = self._args.fasta_guide
        if self._args.bed_guide:
            _check_file(
                self._args.bed_guide,
                self._parser,
                f"Cannot find input guide BED {self._args.bed_guide}",
            )
            self._guidebed = self._args.bed_guide

    def _validate_pam(self) -> None:
        if any(nt.upper() not in IUPAC for nt in self._args.pam):
            self._parser.error(f"Invalid PAM sequence {self._args.pam}")

    def _validate_mm(self) -> None:
        if self._args.mm < 0:
            self._parser.error(
                f"Invalid number of mismatches selected ({self._args.mm})"
            )

    def _validate_bdna(self) -> None:
        if self._args.bdna < 0:
            self._parser.error(
                f"Invalid number of DNA bulges selected ({self._args.bdna})"
            )

    def _validate_brna(self) -> None:
        if self._args.brna < 0:
            self._parser.error(
                f"Invalid number of RNA bulges selected ({self._args.brna})"
            )

    def _validate_annotation(self) -> None:
        for bed in self._args.annotations:
            _check_file(bed, self._parser, f"Cannot find annotation file {bed}")
            if not bed.endswith((".bed", ".bed.gz")):
                self._parser.error(
                    f"Unsupported annotation file '{bed}'. "
                    "Expected a .bed or .bed.gz file"
                )
        self._annotations: List[str] = self._args.annotations
        if len(self._annotations) > 10:
            self._parser.error(
                f"Too many input annotation files: {len(self._annotations)}. "
                "Maximum number of supported by {} annotation files is 10"
            )

    def _validate_annotation_names(self) -> None:
        names: List[str] = self._args.annotation_colnames
        if names is not None and len(names) != len(self._args.annotations):
            self._parser.error(
                f"Number of --annotation-names ({len(names)}) does not match "
                f"the number of --annotations files ({len(self._args.annotations)})"
            )
        self._annotation_names = (
            list(names)
            if names
            else [f"annotation_{i}" for i, _ in enumerate(names, start=1)]
        )

    def _check_consistency(self) -> None:
        """Check the consistency and validity of parsed input arguments.

        Validates the existence, type, and content of input files and directories,
        and sets the list of VCF files found in the specified directory.

        Returns:
            None
        """
        self._validate_genome_folder()
        if self._args.vcf:
            self._validate_vcf_folder()
        self._validate_guides()
        self._validate_pam()
        self._validate_mm()
        self._validate_bdna()
        self._validate_brna()
        if self._args.annotations:
            self._validate_annotation()
            self._validate_annotation_names()
        self._validate_output_folder()
        self._validate_threads()

    @property
    def fastas(self) -> List[str]:
        return self._fastas

    @property
    def vcfs(self) -> List[str]:
        return self._vcfs if hasattr(self, "_vcfs") else []

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
        return self._args.mm

    @property
    def bdna(self) -> int:
        return self._args.bdna

    @property
    def brna(self) -> int:
        return self._args.brna

    @property
    def annotations(self) -> List[str]:
        return self._annotations if hasattr(self, "_annotations") else []

    @property
    def annotation_names(self) -> List[str]:
        return self._annotation_names if hasattr(self, "_annotation_names") else []

    @property
    def upstream(self) -> bool:
        return self._args.upstream

    @property
    def outdir(self) -> str:
        return self._outdir

    @property
    def threads(self) -> int:
        return self._threads


# ==============================================================================
# Internal helpers
# ==============================================================================


def _check_folder(dirname: str, parser: Crisprme2ArgumentParser, msg: str) -> None:
    """Report a usage error if *dirname* is not an existing directory.

    Parameters
    ----------
    dirname : str
        Directory path to check.
    parser : Crisprme2ArgumentParser
        Parser used to report the error.
    msg : str
        Error message shown on failure.

    Returns
    -------
    None
    """
    if not os.path.exists(dirname) or not os.path.isdir(dirname):
        parser.error(msg)


def _check_file(fname: str, parser: Crisprme2ArgumentParser, msg: str) -> None:
    """Report a usage error if *fname* is not an existing file.

    Parameters
    ----------
    fname : str
        File path to check.
    parser : Crisprme2ArgumentParser
        Parser used to report the error.
    msg : str
        Error message shown on failure.

    Returns
    -------
    None
    """
    if not os.path.exists(fname) or not os.path.isfile(fname):
        parser.error(msg)


def _validate_threads_num(threads: int, parser: Crisprme2ArgumentParser) -> int:
    """Validate a thread count against the available CPU cores.

    A value of ``0`` is interpreted as "use all cores".

    Parameters
    ----------
    threads : int
        Requested number of threads.
    parser : Crisprme2ArgumentParser
        Parser used to report the error.

    Returns
    -------
    int
        The validated thread count (all cores when *threads* is ``0``).

    Raises
    ------
    SystemExit
        Via the parser's ``error`` method if *threads* is negative or exceeds
        the available cores.
    """
    max_threads = multiprocessing.cpu_count()
    if threads < 0 or threads > max_threads:
        parser.error(
            f"Forbidden number of threads provided ({threads}). "
            f"Max number of available cores: {max_threads}"
        )
    return max_threads if threads == 0 else threads


def _check_retrieved_files(
    fnames: List[str], parser: Crisprme2ArgumentParser, msg: str
) -> None:
    """Report a usage error if *fnames* is empty.

    Parameters
    ----------
    fnames : List[str]
        Collected file paths.
    parser : Crisprme2ArgumentParser
        Parser used to report the error.
    msg : str
        Error message shown on failure.

    Returns
    -------
    None
    """
    if not fnames:
        parser.error(msg)
