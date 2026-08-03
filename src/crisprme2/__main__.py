"""
CRISPRme2 {version}

Copyright (C) 2025 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>

CRISPRme2: High-performance and scalable tool for variant- and haplotype-aware genome-wide off-target
assessment in CRISPR-Cas systems

CRISPRme2 is a high-performance and scalable tool for genome-wide off-target assessment in CRISPR-Cas
systems. It supports variant-aware and haplotype-aware predictions, integrating SNVs, indels, and
population-scale haplotypes with orthogonal genomic annotations to prioritize off-targets across personal
and population genomes


Usage:
    crisprme2 complete-search --genome <genome-dir> --vcf <vcf-dir> --guide <guide>

Run 'crisprme2 -h/--help' to display the complete help
"""

from argparse import _SubParsersAction
from time import time

import sys
import os


from .crisprme2_argparse import Crisprme2ArgumentParser
from .crisprme2_inputargs import Crisprme2SearchInputArgs
from .search import execute_offtargets_search
from .exception_handlers import sigint_handler
from .utils import COMMANDS, TOOLNAME
from .version import __version__


def create_parser_crisprme2() -> Crisprme2ArgumentParser:
    """Creates and configures the main argument parser for the CRISPRme2 CLI.

    This function sets up the command-line interface, including all available
    commands and their arguments, for the CRISPRme2 toolkit.

    Returns:
        Crisprme2ArgumentParser: The configured argument parser for CRISPRme2.
    """
    # force displaying docstring at each usage display and force
    # the default help to not being shown
    parser = Crisprme2ArgumentParser(usage=__doc__, add_help=False)  # type: ignore
    group = parser.add_argument_group("Options")  # arguments group
    # add help and version arguments
    group.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )
    group.add_argument(
        "--version",
        action="version",
        help=f"Show {TOOLNAME} version and exit",
        version=__version__,
    )
    # create subparsers for different functionalities
    subparsers = parser.add_subparsers(
        dest="command",
        title="Available commands",
        metavar="",  # needed for help formatting (avoid <command to be displayed>)
        description=None,
    )
    # crisprme2 complete-search command
    create_search_parser(subparsers)
    return parser


def create_search_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    """Creates the argument parser for the CRISPRme2 complete-search command.

    This function defines and configures all arguments and options available for
    the search functionality of CRISPRme2.

    Args:
        subparser (_SubParsersAction): The subparsers object to which the search
            parser will be added.

    Returns:
        _SubParsersAction: The configured search command parser.
    """
    parser_search = subparser.add_parser(
        COMMANDS[0],
        usage="CRISPRme2 complete-search {version}\n\nUsage:\n"
        "\ncrisprme2 complete-search --genome <genome-dir> --vcf <vcf-dir> "
        "--guide <guide> --pam <pam> --outdir <output-dir>\n\n",
        description="Automated end-to-end search pipeline that processes raw input "
        "data through off-targets identification, scoring, and annotation of results",
        help="perform a comprehensive off-targets search across the reference genome "
        "and optionally variant-aware genomes. Includes CFD, CRISTA (for Cas9 "
        "systems), CRISPR-bulge, and Elevation score (for compatible Cas systems) "
        "to evaluate genetic diversity impact on off-targets, and automated "
        "targets annotation",
        add_help=False,
    )
    general_group = parser_search.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    required_group = parser_search.add_argument_group("Options")
    required_group.add_argument(
        "--genome",
        type=str,
        metavar="GENOME-DIR",
        required=True,
        dest="genome",
        help="folder containing genome FASTA files for off-targets search. Each "
        "chromosome must be in a separate FASTA file (e.g., chr1.fa, chr2.fa). "
        "All files in the folder will be used as the reference genome",
    )
    required_group.add_argument(
        "--pam",
        type=str,
        metavar="PAM",
        required=True,
        dest="pam",
        help="PAM sequence (e.g., NGG, NRG, TTTV, etc.)",
    )
    guide_group = required_group.add_mutually_exclusive_group(required=True)
    guide_group.add_argument(
        "--guide",
        type=str,
        dest="guide",
        metavar="GUIDE",
        help="guide RNA sequence (spacer only, without PAM) used to search for "
        "potential off-targets in both the reference and alternative genomes. "
        "Cannot be used with --sequence or --coordinates",
    )
    guide_group.add_argument(
        "--sequence",
        type=str,
        dest="fasta_guide",
        metavar="FASTA-FILE",
        help="FASTA file containing guide sequences. Cannot be used with --guide "
        "or --coordinates",
    )
    guide_group.add_argument(
        "--coordinates",
        type=str,
        dest="bed_guide",
        metavar="BED-FILE",
        help="BED file with genomic coordinates for guide regions. Cannot be "
        "used with --guide or --sequence",
    )
    required_group.add_argument(
        "--mm",
        type=int,
        metavar="MISMATCHES",
        dest="mm",
        required=True,
        help="maximum number of mismatches allowed between the guide and off-targets",
    )
    required_group.add_argument(
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        nargs="?",
        default=os.getcwd(),
        help="output directory where reports and results will be saved. "
        "(default: current working directory)",
    )
    optional_group = parser_search.add_argument_group("Optional arguments")
    optional_group.add_argument(
        "--vcf",
        type=str,
        metavar="VCF-DIR",
        dest="vcf",
        nargs="?",
        default="",
        help="optional folder storing VCF files to consider in the off-targets search. "
        "(default: no variant-aware analysis)",
    )
    optional_group.add_argument(
        "--bdna",
        type=int,
        dest="bdna",
        metavar="NUM-BULGE-DNA",
        required=False,
        default=0,
        help="maximum number of DNA bulges allowed in the search (default: 0)",
    )
    optional_group.add_argument(
        "--brna",
        type=int,
        dest="brna",
        metavar="NUM-BULGE-RNA",
        required=False,
        default=0,
        help="maximum number of RNA bulges allowed in the search (default: 0)",
    )
    optional_group.add_argument(
        "--max-edit-distance",
        type=int,
        dest="max_edit_dist",
        metavar="MAX-EDIT-DISTANCE",
        required=False,
        default=0,
        help="maximum allowed edit distance between the guide RNA and aligned "
        "target. The edit distance is computed as mismatches + RNA bulges + "
        "DNA bulges. A value of 0 disables this filter and automatically uses "
        "the maximum search distance defined by the mismatch and bulge "
        "thresholds (default: mismatches + RNA bulges + DNA bulges)",
    )
    optional_group.add_argument(
        "--annotation",
        type=str,
        metavar="ANNOTATION-BED",
        dest="annotations",
        nargs="*",
        default=[],
        help="one or more BED files specifying genomic regions used to annotate "
        "guide candidates. Each file should follow the standard BED format "
        "(at least: chrom, start, end), and should include additional annotation "
        "on the 4th column (default: no annotation)",
    )
    optional_group.add_argument(
        "--annotation-colnames",
        type=str,
        metavar="ANNOTATION-COLNAMES",
        dest="annotation_colnames",
        nargs="*",
        default=[],
        help="list of custom column names to use in the final report. Each name "
        "corresponds to one of the input BED files provided with '--annotation'. "
        "Must match the number and order of files in '--annotation' (default: "
        "annotation columns are named 'annotation_<i>')",
    )
    optional_group.add_argument(
        "--upstream",
        action="store_true",
        dest="upstream",
        default=False,
        help="if set, PAM occurs upstream (left side) of the guide "
        "(default: PAM occurs downstream (right side))",
    )
    optional_group.add_argument(
        "--cluster-distance",
        type=int,
        dest="cluster_dist",
        metavar="CLUSTER-DIST",
        required=False,
        default=3,
        help="maximum genomic distance (in bp) between off-target alignments to "
        "group them into a single editing site. Alignments within this distance "
        "are reported as alternative alignments of the same site rather than as "
        "independent off-targets (default: 3)",
    )
    optional_group.add_argument(
        "--prioritization-criteria",
        type=str,
        dest="prioritization_criteria",
        metavar="PRIORITIZATION-CRITERIA",
        required=False,
        default="edit-dist,bdna,brna,mm",
        help="comma-separated list of criteria, in descending priority order, "
        "used to select the primary alignment when multiple alignments belong "
        "to the same editing site. Remaining alignments are reported as "
        "alternative alignments. Available criteria: 'edit-dist' (total edit "
        "distance), 'bdna' (DNA bulges), 'brna' (RNA bulges), 'mm' (mismatches), "
        "'cfd' (CFD score), 'crista' (CRISTA score), 'elevation' (Elevation "
        "score) (default: edit-dist,bdna,brna,mm)",
    )
    optional_group.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        dest="threads",
        nargs="?",
        default=1,
        help="number of threads. Use 0 for using all available cores (default: 1)",
    )
    return parser_search


def main():
    start = time()  # track elapsed time
    try:
        parser = create_parser_crisprme2()  # parse input argument using custom parser
        if not sys.argv[1:]:  # no input args -> print help and exit
            parser.error_noargs()
        args = parser.parse_args(sys.argv[1:])  # parse input args
        if args.command == COMMANDS[0]:  # complete-search command
            execute_offtargets_search(Crisprme2SearchInputArgs(args, parser))
    except KeyboardInterrupt:
        sigint_handler()  # catch SIGINT and exit gracefully
    sys.stdout.write(f"{TOOLNAME} - Elapsed time {(time() - start):.2f}s\n")


# --------------------------------> ENTRY POINT <--------------------------------
if __name__ == "__main__":
    main()
