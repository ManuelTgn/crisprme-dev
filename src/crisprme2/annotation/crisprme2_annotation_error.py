"""
crisprme2_scores_error.py
---------------------
Custom exception hierarchy for the crisprme2 annotation API.

All public exceptions derive from Crisprme2Error so callers can catch the
entire family with a single ``except Crisprme2Error`` clause while still
being able to discriminate between subsystems when needed.

Hierarchy
~~~~~~~~~
::

    Crisprme2Error
    └── Crisprme2AnnotationError                    - Generic annotation error
        └── Crisprme2FunctionalAnnotationError      - Functional annotation
"""

from ..crisprme2_error import Crisprme2Error


class Crisprme2AnnotationError(Crisprme2Error):
    """Raised when the scoring calculation encounters an unrecoverable error"""

    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2FunctionalAnnotationError(Crisprme2AnnotationError):
    """Raised when the CFD scoring calculation encounters an unrecoverable error"""

    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
