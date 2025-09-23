""" """

from .crisprme2_argparse import Crisprme2SearchInputArgs
from .logger import CrisprmeLoggers
from .guide import read_guides
from .utils import TOOLNAME

import os

def complete_search(args: Crisprme2SearchInputArgs) -> None:
    loggers = CrisprmeLoggers()  # initialize loggers
    guides = read_guides(args.guide, "", "", loggers)
    print(guides)
