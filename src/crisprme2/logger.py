""" """

from typing import Optional, Union
from logging.handlers import RotatingFileHandler

import logging
import os

LOGLEVELS = [0, 1, 2]

class BaseLogger:
    def __init__(self, name: str, log_file: str, log_dir: str = "logs", level: int = logging.INFO, max_bytes : int = 10**6, backup_count: int = 5) -> None:
        self.logger = logging.getLogger(name + "_" + log_file)
        self.logger.setLevel(logging.DEBUG)  # Capture all; filter later

        if not self.logger.handlers:
            os.makedirs(log_dir, exist_ok=True)
            file_path = os.path.join(log_dir, log_file)
            handler = RotatingFileHandler(
                file_path, maxBytes=max_bytes, backupCount=backup_count  
            )
            handler.setLevel(level)  
            handler.setFormatter(self._get_formatter())
            self.logger.addHandler(handler)
            log_filter = self._get_filter(level)  
            if log_filter:
                handler.addFilter(log_filter)

    def _get_formatter(self) -> logging.Formatter:
        return logging.Formatter(
            fmt="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )

    def _get_filter(self, level: int) -> Optional[logging.Filter]:  # override in subclasses if needed
        return None

    def get_logger(self) -> logging.Logger:
        return self.logger


class BasicLogger(BaseLogger):
    def __init__(self, name: str, log_dir: str = "logs"):
        super().__init__(name, log_file="basic.log", log_dir=log_dir, level=logging.INFO)

    def _get_filter(self, level: int):
        class InfoOnlyFilter(logging.Filter):
            def filter(self, record):
                return record.levelno == logging.INFO
        return InfoOnlyFilter()


class VerboseLogger(BaseLogger):
    def __init__(self, name: str, log_dir: str = "logs"):
        super().__init__(name, log_file="verbose.log", log_dir=log_dir, level=logging.DEBUG)

    def _get_filter(self, level: int):
        class DebugAndInfoFilter(logging.Filter):
            def filter(self, record):
                return record.levelno in (logging.DEBUG, logging.INFO)
        return DebugAndInfoFilter()


class ErrorLogger(BaseLogger):
    def __init__(self, name: str, log_dir: str = "logs"):
        super().__init__(name, log_file="errors.log", log_dir=log_dir, level=logging.ERROR)


def logger(log_file: logging.Logger, message: str, level: int) -> None:
    assert level in LOGLEVELS  # shouldn't be different
    if level == LOGLEVELS[0]:  # basic log
        log_file.info(message)
    elif level == LOGLEVELS[1]:  # verbose log
        log_file.debug(message)
    else:  # error log
        log_file.error(message)
