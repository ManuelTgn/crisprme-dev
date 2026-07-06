"""
complete_search.py
------------------
Composition root for the CRISPRme2 complete-search pipeline.

Responsibilities
~~~~~~~~~~~~~~~~
This module is the **only** place in the codebase that knows about both
the CLI argument namespace and the internal pipeline components.  It:

1. Constructs domain objects from raw CLI values (PAM, guides, thresholds).
2. Assembles the ordered list of transform callables (scorers, annotators).
3. Calls the search entry-point with fully-wired arguments.

What this module does NOT do
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- It does not parse CLI arguments — that is ``__main__``'s responsibility.
- It does not implement pipeline logic — that is ``search.py``'s responsibility.
- It does not implement scoring — that is each scorer's responsibility.

Adding a new transform
~~~~~~~~~~~~~~~~~~~~~~
To add a new scoring or annotation transform to the pipeline:

1. Implement the :class:`~crisprme2.protocol.Transformer` protocol in the
   appropriate subpackage (e.g. ``crisprme2.scores.crista``).
2. Instantiate it in :func:`_build_transforms` and append it to the list.
3. Nothing else needs to change.
"""

from __future__ import annotations

from typing import List, Tuple

from .crisprme_core_api import Thresholds
from .crisprme2_argparse import Crisprme2SearchInputArgs
from .crisprme2 import TOOLNAME
from .guide import read_guides, GuidesList
from .logger import CrisprmeLoggers
from .pam import read_pam, PAM
from .protocol import Transformer
from .scores import CfdScorer
from .search import search_offtargets_reference_genome


def _build_pam_and_guides(
    args: Crisprme2SearchInputArgs, loggers: CrisprmeLoggers
) -> Tuple[GuidesList, PAM]:
    """
    Initialise PAM and guide data structures from validated CLI arguments.

    Parameters
    ----------
    args : Crisprme2SearchInputArgs
        Validated argument namespace.
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Returns
    -------
    tuple[GuidesList, PAM]
        ``(guides, pam)`` ready for use in the search pipeline.
    """
    loggers.basiclog.info("Initialising PAM and guide data structures")
    pam = read_pam(args.pam, loggers)
    guides = read_guides(args, loggers)
    loggers.verboselog.debug(f"PAM: {pam} | guides: {len(guides)}")
    return guides, pam


def _build_thresholds(
    args: Crisprme2SearchInputArgs, loggers: CrisprmeLoggers
) -> Thresholds:
    """
    Construct a :class:`~crisprme2.crisprme_core_api.Thresholds` instance
    from validated CLI arguments.

    Parameters
    ----------
    args : Crisprme2SearchInputArgs
        Validated argument namespace.
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Returns
    -------
    Thresholds
        Alignment thresholds for this run.
    """
    loggers.verboselog.debug(
        f"Building Thresholds(max_mm={args.mm}, bdna={args.bdna}, brna={args.brna})"
    )
    return Thresholds(
        max_mm=args.mm, max_bdna=args.bdna, max_brna=args.brna, loggers=loggers
    )


class ExampleTransform:
    def __call__(self):
        pass


def _build_transforms(pam: PAM, loggers: CrisprmeLoggers) -> List[Transformer]:
    transforms: List[Transformer] = []
    # ---- scoring transform
    # CFD score + slot 0
    # CFD pam is the last two bases of the PAM sequence
    # For NGG the key is "GG"; for NGA it is "GA", etc.
    pam_key = pam.pam[-2:]
    transforms.append(CfdScorer(pam=pam_key, loggers=loggers))

    # ---> future scorers <---

    loggers.verboselog.debug(
        "Transform chain assembled: " f"{[type(t).__name__ for t in transforms]}"
    )
    return transforms


def execute_complete_search(args: Crisprme2SearchInputArgs) -> None:
    """
    Run the full CRISPRme2 complete-search pipeline.

    This is the composition root: it wires CLI arguments to pipeline
    components and delegates execution to specialised modules.  The call
    graph is::

        execute_complete_search(args)
            ├── CrisprmeLoggers(args.outdir)
            ├── _build_pam_and_guides(args)     -> GuidesList, PAM
            ├── _build_thresholds(args)         -> Thresholds
            ├── _build_transforms(pam)          -> list[Transformer]
            └── (per guide)
                └── search_offtargets_reference_genome(...)

    Parameters
    ----------
    args : Crisprme2SearchInputArgs
        Fully validated CLI argument namespace produced by
        :func:`~crisprme2.__main__.create_parser_crisprme2`.

    Raises
    ------
    Crisprme2SearchError
        If any component of the search pipeline fails.
    """
    loggers = CrisprmeLoggers(args.outdir)  # initialize loggers
    loggers.basiclog.info(f"Start {TOOLNAME} search")

    # initialize pam and guide objects
    guides, pam = _build_pam_and_guides(args, loggers)
    # initialize thresholds object
    thresholds = _build_thresholds(args, loggers)
    # initialize transforms
    transforms = _build_transforms(pam, loggers)

    for guide in guides:
        # retrieve candidate off-targets for current guide
        loggers.verboselog.debug(
            f"Starting off-target search for guide {guide.sequence}"
        )
        if args.vcfs:
            # variant and haplotype aware search path (not yet implemented)
            loggers.verboselog.debug(
                "VCF files provided - variant-aware search path "
                "not yet implemented (skipping)"
            )
            continue
        search_offtargets_reference_genome(
            args.fastas,
            pam,
            guide,
            args.upstream,
            args.threads,
            thresholds,
            transforms,
            loggers,
        )
