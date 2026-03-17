""" """

from ..crisprme2_error import Crisprme2Error


class Crisprme2BatcherError(Crisprme2Error):
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


class Crisprme2AlignmentError(Crisprme2Error):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
