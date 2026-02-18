""" """

from crisprme2.logger import CrisprmeLoggers
from crisprme2.crisprme2_error import Crisprme2CfdScoreError
from crisprme2.sequence import dna2rna, reverse_complement

from typing import Dict, Tuple

import pickle
import os


# define cfd models folder name
MODELSDIR = "models"

# define mismatch scores pickle file name
MMSCORES = "mismatch_score.pkl"

# define pam scores pickle file name
PAMSCORES = "pam_scores.pkl"


def load_mismatch_pam_scores(loggers: CrisprmeLoggers) -> Tuple[Dict[str, float], Dict[str, float]]:
    modelspath = os.path.join(os.path.abspath(os.path.dirname(__file__)), MODELSDIR)
    try:  # load mismatches and PAM scores (Doench et al., 2016)
        mmscores = pickle.load(open(os.path.join(modelspath, MMSCORES), mode="rb"))
        pamscores = pickle.load(open(os.path.join(modelspath, PAMSCORES), mode="rb"))
    except (pickle.UnpicklingError, pickle.PickleError, Exception) as e:
        loggers.errorlog.log_raise_exception(f"Failed loading CFD models: {e}", os.EX_IOERR, Crisprme2CfdScoreError)
    return mmscores, pamscores

def _construct_cfd_key(wt_nt: str, nt_sg: str, pos: int, loggers: CrisprmeLoggers) -> str:
    nt_sg = reverse_complement(nt_sg, loggers)
    return f"r{wt_nt}:d{nt_sg},{pos + 1}"


def _compute_cfd_mismatches(wildtype: str, sg: str, mmscores: Dict[str, float], loggers: CrisprmeLoggers) -> float:
    try:# weight cfd score by mismatches using cfd model
        score = 1.0  # initialize cfd score
        for i, nt_sg in enumerate(sg):
            if i >= 20:  # handle presence of bulges in off-targets
                break
            if wildtype[i] == nt_sg:  # guide-target match
                score *= 1
                continue
            elif wildtype[i] == "-" or nt_sg == "-":  # bulge encountered
                score *= 1
                continue
            # query cfd mismatch scoring model 
            keymap = _construct_cfd_key(wildtype[i], nt_sg, i, loggers)
            score *= mmscores[keymap]
    except (KeyError, ValueError) as e:
        loggers.errorlog.log_raise_exception(f"Failed CFD mismatch weighting: {e}", os.EX_DATAERR, Crisprme2CfdScoreError)
    return score

def _compute_cfd_pam(score: float, pam: str, pamscores: Dict[str, float], loggers: CrisprmeLoggers) -> float:
    try:  # weight cfd by pam score
        return score * pamscores[pam]  
    except (KeyError, ValueError) as e:
        loggers.errorlog.log_raise_exception(f"Failed CFD PAM weighting: {e}", os.EX_DATAERR, Crisprme2CfdScoreError)


def compute_cfd(wildtype: str, sg: str, pam: str, mmscores: Dict[str, float], pamscores: Dict[str, float], loggers: CrisprmeLoggers) -> float:
    # convert off-target and guide sequence to rna
    wildtype, sg = dna2rna(wildtype, loggers).upper(), dna2rna(sg, loggers).upper()
    # compute cfd score using mismatch weights
    score = _compute_cfd_mismatches(wildtype, sg, mmscores, loggers)
    # weight cfd score by pam score
    score = _compute_cfd_pam(score, pam.upper(), pamscores, loggers)
    return score

        



