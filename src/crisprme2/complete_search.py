""" """

from .crisprme2_argparse import Crisprme2SearchInputArgs
from .crisprme2 import TOOLNAME
from .logger import CrisprmeLoggers
from .search import search_offtargets_reference_genome
from .guide import read_guides, GuidesList, Guide
from .pam import read_pam, PAM

from typing import Tuple, List


def _init_pam_guide(
    args: Crisprme2SearchInputArgs, loggers: CrisprmeLoggers
) -> Tuple[GuidesList, PAM]:
    loggers.basiclog.info("Initialize PAM and guides data structures")
    pam = read_pam(args.pam, loggers)  # pam object
    guides = read_guides(args, loggers)  # guides list object
    return guides, pam


def _retrieve_target_candidates(
    fastas: List[str],
    vcfs: List[str],
    pam: PAM,
    guide: Guide,
    offset: int,
    right: bool,
    threads: int,
    loggers: CrisprmeLoggers,
):
    if vcfs:  # variant- and haplotype-aware search
        pass  # insert here call to enrichment module
    else:  # reference/assembly fasta search
        loggers.verboselog.debug(
            "Reference/assembly genome off-targets extraction pipeline"
        )
        search_offtargets_reference_genome(fastas, pam, guide, offset, right, threads, loggers)


def execute_complete_search(args: Crisprme2SearchInputArgs) -> None:

    # engine = init_engine(args)

    loggers = CrisprmeLoggers(args.outdir)  # initialize crisprme loggers
    loggers.basiclog.info(f"Start {TOOLNAME} search")
    # initialize guides and pam objects
    guides, pam = _init_pam_guide(args, loggers)
    # retrieve candidate off-targets for each guide
    loggers.basiclog.info("Retrieve candidate off-targets")
    for guide in guides:
        loggers.verboselog.debug(f"Retrieve off-targets for guide {guide}")
        _retrieve_target_candidates(
            args.fastas,
            args.vcfs,
            pam,
            guide,
            max(args.bdna, args.brna),
            args.right,
            args.threads,
            loggers,
            # engine
        )
