""" """

from .utils import TOOLNAME, warning

from typing import Optional, Dict, NoReturn
from logging.handlers import RotatingFileHandler
from colorama import Fore, Style

import logging
import shutil
import sys
import os


# define log levels
LOGLEVELS = [0, 1, 2]

# define logs directory
LOGDIR = "logs"


class BaseLogger:

    _logs_cleaned = False  # class variable to track if logs have been cleaned this run

    def __init__(
        self,
        name: str,
        log_file: str,
        log_dir: str = "logs",
        level: int = logging.INFO,
        max_bytes: int = 10**6,
        backup_count: int = 5,
        clean_logs: bool = True,
    ) -> None:
        # clean logs folder only once per run
        if clean_logs and not BaseLogger._logs_cleaned:
            self._clean_logs_folder(log_dir)
            BaseLogger._logs_cleaned = True  # logs folder is clean now
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

    def _clean_logs_folder(self, log_dir: str) -> None:
        if os.path.exists(log_dir) and os.path.isdir(log_dir):
            try:
                for fname in os.listdir(
                    log_dir
                ):  # remove all files in the logs folders
                    fpath = os.path.join(log_dir, fname)
                    try:
                        if os.path.isfile(fpath) or os.path.islink(fpath):
                            os.unlink(fpath)  # delete log
                        elif os.path.isdir(fpath):  # what? better remove it anyway
                            shutil.rmtree(fpath)
                    except Exception as e:
                        warning(f"Failed to delete {fpath}. Reason {e}", 1)
            except Exception as e:  # throw a warning, but do not halt execution
                warning(f"Failed to clean logs folder {log_dir}. Reason: {e}", 1)

    @classmethod
    def reset_cleanup_flag(cls):
        cls._logs_cleaned = False

    def _get_formatter(self) -> logging.Formatter:
        return logging.Formatter(
            fmt="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    def _get_filter(
        self, level: int
    ) -> Optional[logging.Filter]:  # override in subclasses if needed
        return None

    def get_logger(self) -> logging.Logger:
        return self.logger


class BasicLogger(BaseLogger):
    def __init__(self, name: str, log_dir: str = "logs"):
        super().__init__(
            name, log_file="basic.log", log_dir=log_dir, level=logging.INFO
        )

    def _get_filter(self, level: int):
        class InfoOnlyFilter(logging.Filter):
            def filter(self, record):
                return record.levelno == logging.INFO

        return InfoOnlyFilter()

    def info(self, message: str) -> None:
        self.logger.info(message)  # write basic info level to log


class VerboseLogger(BaseLogger):
    def __init__(self, name: str, log_dir: str = "logs"):
        super().__init__(
            name, log_file="verbose.log", log_dir=log_dir, level=logging.DEBUG
        )

    def _get_filter(self, level: int):
        class DebugAndInfoFilter(logging.Filter):
            def filter(self, record):
                return record.levelno in (logging.DEBUG, logging.INFO)

        return DebugAndInfoFilter()

    def info(self, message: str) -> None:
        self.logger.info(message)  # write basic info level to log

    def debug(self, message: str) -> None:
        self.logger.debug(message)  # write debug level to log


class ErrorLogger(BaseLogger):
    def __init__(self, name: str, log_dir: str = "logs"):
        super().__init__(
            name, log_file="errors.log", log_dir=log_dir, level=logging.ERROR
        )
        self._log_dir = log_dir

    def log_exception(self, message: str, code: int, exc_info: bool = True) -> NoReturn:
        _halt_message(self._log_dir)  # print execution halt message in terminal
        self.logger.error(message, exc_info=exc_info)
        sys.exit(code)  # halt execution

    def log_raise_exception(
        self,
        message: str,
        code: int,
        exception_type: type = Exception,
        exc_info: bool = True,
    ) -> NoReturn:
        try:  # force exception raise to capture it in log file
            raise exception_type(message)
        except exception_type:  # type: ignore
            self.log_exception(message, code, exc_info=exc_info)

    def log_error_with_context(
        self, message: str, code: int, context: Optional[Dict] = None
    ) -> NoReturn:
        if context:
            context_str = " | ".join([f"{k}={v}" for k, v in context.items()])
            fullmessage = f"{message} | Context: {context_str}"
        else:
            fullmessage = message
        _halt_message(self._log_dir)  # print execution halt message in terminal
        self.logger.error(fullmessage)
        sys.exit(code)  # halt execution


class CrisprmeLoggers:
    def __init__(self, outdir: str) -> None:
        log_dir = os.path.join(outdir, LOGDIR)
        self._basiclog = BasicLogger(TOOLNAME, log_dir=log_dir)  # 1) basic run info
        self._verboselog = VerboseLogger(
            TOOLNAME, log_dir=log_dir
        )  # 2) verbose debug + info
        self._errorlog = ErrorLogger(TOOLNAME, log_dir=log_dir)  # 3) error + critical

    @property
    def basiclog(self) -> BasicLogger:
        return self._basiclog

    @property
    def verboselog(self) -> VerboseLogger:
        return self._verboselog

    @property
    def errorlog(self) -> ErrorLogger:
        return self._errorlog


def _halt_message(log_dir: str) -> None:
    log_dir = os.path.abspath(log_dir)  # absolute path to logs folder
    haltmsg = (
        f"{Fore.RED}{Style.BRIGHT}\nERROR: {TOOLNAME} run failed.\n\n"
        f"Check log files in: {log_dir} for more information.{Fore.RESET}\n\n"
    )
    sys.stderr.write(haltmsg)
