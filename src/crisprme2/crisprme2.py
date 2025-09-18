""" """

from .crisprme2_argparse import Crisprme2SearchInputArgs
from .logger import BasicLogger, VerboseLogger, ErrorLogger, logger
from .utils import TOOLNAME

def complete_search(args: Crisprme2SearchInputArgs) -> None:
    # initialize loggers:
    basiclog = BasicLogger(TOOLNAME)  #  1) basic info
    verboselog = VerboseLogger(TOOLNAME)  # 2) verbose debug+info
    errorlog = ErrorLogger(TOOLNAME)  # 3) error error+critical
    print("starting search...")