""" """

from abc import ABC, abstractmethod
from argparse import Namespace
from glob import glob
from typing import List, Optional

import multiprocessing
import os


from .assembly import AssemblyInputs
from .crisprme2_argparse import Crisprme2ArgumentParser
from .crisprme2_error import Crisprme2AssemblyError
from .dna_alphabet import DNA, IUPAC
from .utils import CRITERIA


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
        threads: int = self._args.threads
        max_threads = multiprocessing.cpu_count()
        if threads < 0 or threads > max_threads:
            self._parser.error(
                f"Forbidden number of threads provided ({threads}). "
                f"Max number of available cores: {max_threads}"
            )
        self._threads = max_threads if threads == 0 else threads

    @property
    def outdir(self) -> str:
        return self._outdir

    @property
    def threads(self) -> int:
        return self._threads


class Crisprme2SearchInputArgsBase(Crisprme2InputArgs, ABC):

    def __init__(self, args: Namespace, parser: Crisprme2ArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()

    @abstractmethod
    def _validate_inputs(self) -> None:
        raise NotImplementedError

    def _check_consistency(self) -> None:
        self._validate_inputs()  # genome/vcf OR assemblies/chains
        self._validate_guides()
        self._validate_pam()
        self._validate_mm()
        self._validate_bdna()
        self._validate_brna()
        self._validate_max_edit_dist()
        if self._args.annotations:
            self._validate_annotation()
            self._validate_annotation_names()
        self._validate_cluster_dist()
        self._validate_prioritization_criteria()
        self._validate_output_folder()
        self._validate_output_prefix()
        self._validate_threads()

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

    def _validate_max_edit_dist(self) -> None:
        max_edit_dist: int = self._args.max_edit_dist
        if max_edit_dist < 0:
            self._parser.error(
                f"Negative maximum edit distance value selected: {max_edit_dist}"
            )
        max_edit_dist_fx = self.mm + self.bdna + self.brna
        if max_edit_dist > max_edit_dist_fx:
            self._parser.error(
                f"Maximum edit distance value ({max_edit_dist}) above maximum threshold ({max_edit_dist_fx})"
            )
        self._max_edit_dist = max_edit_dist

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

    def _validate_cluster_dist(self) -> None:
        cluster_dist: int = self._args.cluster_dist
        if cluster_dist <= 0:
            self._parser.error(
                f"Cluster distance cannot be 0 or negative: {cluster_dist}"
            )
        self._cluster_dist = cluster_dist

    def _validate_prioritization_criteria(self) -> None:
        prioritization_criteria: List[str] = list(
            self._args.prioritization_criteria.split(",")
        )
        if any(c not in CRITERIA for c in prioritization_criteria):
            forbidden_criteria = ",".join(
                c for c in prioritization_criteria if c not in CRITERIA
            )
            self._parser.error(
                f"Forbidden prioritization criteria among input: {forbidden_criteria}"
            )
        self._prioritization_criteria = prioritization_criteria

    def _validate_output_prefix(self) -> None:
        prefix = self._args.output_prefix
        if prefix is None:
            self._output_prefix = None
            return
        prefix = prefix.strip()
        if not prefix:
            self._parser.error("Empty --output-prefix provided")
        if os.sep in prefix or (os.altsep and os.altsep in prefix):
            self._parser.error(
                f"--output-prefix must be a filename prefix, not a path: {prefix!r}"
            )
        self._output_prefix = prefix

    @property
    def guide(self) -> Optional[str]:
        return self._guide

    @property
    def fasta_guide(self) -> Optional[str]:
        return self._guidefasta

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
    def max_edit_dist(self) -> int:
        if self._max_edit_dist == 0:
            return self.mm + self.bdna + self.brna
        return self._max_edit_dist

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
    def cluster_dist(self) -> int:
        return self._cluster_dist

    @property
    def prioritization_criteria(self) -> List[str]:
        return self._prioritization_criteria

    @property
    def output_prefix(self) -> Optional[str]:
        return self._output_prefix


class Crisprme2SearchInputArgs(Crisprme2SearchInputArgsBase):

    def __init__(self, args: Namespace, parser: Crisprme2ArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()

    def _validate_inputs(self) -> None:
        self._validate_genome_folder()
        if self._args.vcf:
            self._validate_vcf_folder()

    def _validate_genome_folder(self) -> None:
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
        _check_folder(
            self._args.vcf, self._parser, f"Cannot find VCF folder {self._args.vcf}"
        )
        self._vcfs = glob(os.path.join(self._args.vcf, "*.vcf.gz"))
        _check_retrieved_files(
            self._vcfs, self._parser, f"No VCF file found in {self._args.vcf}"
        )

    @property
    def fastas(self) -> List[str]:
        return self._fastas

    @property
    def vcfs(self) -> List[str]:
        return self._vcfs if hasattr(self, "_vcfs") else []


class Crisprme2AssemblySearchInputArgs(Crisprme2SearchInputArgsBase):

    def _validate_inputs(self) -> None:
        _check_folder(
            self._args.assemblies,
            self._parser,
            f"Cannot find assemblies folder {self._args.assemblies}",
        )
        _check_folder(
            self._args.chains,
            self._parser,
            f"Cannot find chains folder {self._args.chains}",
        )
        try:
            self._assemblies = AssemblyInputs.discover(
                self._args.assemblies, self._args.chains
            )
        except Crisprme2AssemblyError as e:
            self._parser.error(str(e))

    @property
    def assemblies(self) -> AssemblyInputs:
        return self._assemblies


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
