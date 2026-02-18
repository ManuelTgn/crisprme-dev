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

from .crisprme2_argparse import Crisprme2ArgumentParser, Crisprme2SearchInputArgs
from .complete_search import execute_complete_search
from .exception_handlers import sigint_handler
from .version import __version__
from .crisprme2 import TOOLNAME

from argparse import _SubParsersAction
from time import time

import sys
import os

# crisprhawk commands
SEARCH = "complete-search"
COMMANDS = [SEARCH]


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
        SEARCH,
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
        dest="genome_dir",
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
    guide_group = parser_search.add_mutually_exclusive_group(required=True)
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
        "-o",
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
        "--right",
        action="store_true",
        dest="right",
        default=False,
        help="if set, guides occur downstream (right side) of the PAM "
        "(default: guides occur upstream (left side))",
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
        if args.command == SEARCH:  # complete-search command
            execute_complete_search(Crisprme2SearchInputArgs(args, parser))
    except KeyboardInterrupt:
        sigint_handler()  # catch SIGINT and exit gracefully
    sys.stdout.write(f"{TOOLNAME} - Elapsed time {(time() - start):.2f}s\n")


# --------------------------------> ENTRY POINT <--------------------------------
if __name__ == "__main__":
    main()
