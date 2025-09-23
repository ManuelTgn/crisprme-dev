""" """

from .crisprme2_argparse import Crisprme2SearchInputArgs
from .logger import CrisprmeLoggers
from .guide import read_guides
from .pam import read_pam


def complete_search(args: Crisprme2SearchInputArgs) -> None:
    loggers = CrisprmeLoggers()  # initialize loggers
    # initialize guides and pam objects
    pam = read_pam(args.pam, loggers)
    guides = read_guides(
        args.guide, args.fasta_guide, args.bed_guide, pam, args.right, loggers
    )
    print(guides)
    print(pam)
