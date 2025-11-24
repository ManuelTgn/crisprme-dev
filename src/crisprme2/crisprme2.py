""" """

from .crisprme2_argparse import Crisprme2SearchInputArgs
from .utils import TOOLNAME
from .enricher import retrieve_target_candidates
from .logger import CrisprmeLoggers
from .guide import read_guides, GuidesList
from .pam import read_pam, PAM

from typing import Tuple

def _init_pam_guide(args: Crisprme2SearchInputArgs, loggers: CrisprmeLoggers) -> Tuple[GuidesList, PAM]:
    loggers.basiclog.info("Initialize PAM and guides data structures")
    pam = read_pam(args.pam, loggers)  # pam object
    guides = read_guides(args, loggers)  # guides list object
    return guides, pam


def complete_search(args: Crisprme2SearchInputArgs) -> None:
    loggers = CrisprmeLoggers()  # initialize crisprme loggers
    loggers.basiclog.info(f"Start {TOOLNAME} search")
    # initialize guides and pam objects
    guides, pam = _init_pam_guide(args, loggers)
    # enrich genome with input variants, and retrieve target candidates
    retrieve_target_candidates(args, pam, len(guides[0]), loggers)



    # enrich_genome(args, loggers)
    # assumes all guides share the same length
    # process_genome(args.fastas, pam, len(guides[0]), args.right, args.outdir, loggers)
