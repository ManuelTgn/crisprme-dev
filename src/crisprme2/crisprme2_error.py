"""Custom exception classes for CRISPRme2 error handling.

This module defines a hierarchy of exception classes for specific error types
encountered in the CRISPRme2 tool, enabling precise and descriptive error
reporting throughout the codebase.
"""


class Crisprme2Error(Exception):
    def __init__(self, value: str):
        # initialize exception object when raised
        self._value = value  # error message or error related info

    def __str__(self):
        return repr(self._value)  # string representation for the exception


class Crisprme2FastaError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2GuideError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2PamError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
