"""
crisprme2_scores_error.py
---------------------
Custom exception hierarchy for the crisprme2 scores API.

All public exceptions derive from Crisprme2Error so callers can catch the
entire family with a single ``except Crisprme2Error`` clause while still
being able to discriminate between subsystems when needed.

Hierarchy
~~~~~~~~~
::

    Crisprme2Error
    └── Crisprme2ScoreError                 - Generic score error / runtime failures
            └── Crisprme2CfdScoreError      - CFD score model
"""

from ..crisprme2_error import Crisprme2Error


class Crisprme2ScoreError(Crisprme2Error):
    """Raised when the scoring calculation encounters an unrecoverable error"""

    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2CfdScoreError(Crisprme2ScoreError):
    """Raised when the CFD scoring calculation encounters an unrecoverable error"""

    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
