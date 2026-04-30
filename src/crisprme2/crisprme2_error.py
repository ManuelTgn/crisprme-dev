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


class Crisprme2SequenceError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2ContigSequenceError(Crisprme2SequenceError):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2ReverseComplementError(Crisprme2SequenceError):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2DnaRnaError(Crisprme2SequenceError):
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


class Crisprme2VCFError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2VCFFileNotFoundError(Crisprme2VCFError):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2VCFFormatError(Crisprme2VCFError):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2SampleError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2VariantRecordError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2IupacTableError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2EnrichmentError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2TargetError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2SearchError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2CallbackError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2ScoreError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2CfdScoreError(Crisprme2ScoreError):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2AnnotationError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
