""" """

from .crisprme2_argparse import Crisprme2SearchInputArgs
from .utils import TOOLNAME
from .enricher import enrich_genome
from .logger import CrisprmeLoggers
from .guide import read_guides
from .pam import read_pam


def complete_search(args: Crisprme2SearchInputArgs) -> None:
    loggers = CrisprmeLoggers()  # initialize loggers
    loggers.basiclog.info(f"Start {TOOLNAME} search")
    # initialize guides and pam objects
    pam = read_pam(args.pam, loggers)
    guides = read_guides(args.guide, args.fasta_guide, args.bed_guide, loggers)


    # enrich_genome(args, loggers)
    # assumes all guides share the same length
    # process_genome(args.fastas, pam, len(guides[0]), args.right, args.outdir, loggers)
