""" """

from .crisprme2_error import Crisprme2SampleError
from .logger import CrisprmeLoggers

import os


class Sample:
    def __init__(self, name: str, idx: int, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        if not name or not isinstance(name, str):
            self._loggers.errorlog.log_raise_exception(
                "Sample name must be a non-empty string",
                os.EX_DATAERR,
                Crisprme2SampleError,
            )
        self._name = name  # sample name
        self._idx = idx  # sample identifier

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; name={self._name}"

    @property
    def name(self) -> str:
        return self._name
    
class SampleList:
    pass
