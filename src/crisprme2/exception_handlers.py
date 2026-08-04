"""Exception handling utilities for CRISPRme2.

This module provides the two functions that form the error-handling backbone
of the entire CRISPRme2 pipeline.  Every module that needs to report a
fatal error or respond to a keyboard interrupt delegates to one of these
functions rather than raising exceptions or calling ``sys.exit`` directly,
ensuring that error output is consistent, coloured, and controlled by a single
``debug`` flag.

Design rationale
----------------
CRISPRme2 distinguishes between two runtime modes:

* **Normal mode** (``debug=False``) — errors are presented as concise,
  user-friendly messages written to *stderr* in red (via ``colorama``) and the
  process exits immediately with an appropriate POSIX exit code.  No traceback
  is shown.  This is the default experience for end users running the tool from
  the command line.

* **Debug mode** (``debug=True``) — the typed exception is raised with its
  full message and, when an originating exception is available, exception
  chaining is used (``raise ... from e``) so that the complete traceback is
  preserved for diagnosis.  This mode is activated by the ``--debug`` flag and
  is intended for developers and for automated test suites that ``pytest``
  captures.

Both functions use ``colorama`` for cross-platform ANSI colour support.
:func:`exception_handler` initialises ``colorama`` on every call via
``colorama.init()``; this is lightweight and idempotent, making explicit
package-level initialisation unnecessary.

Typical usage example
---------------------
::

    from crisprme2.exception_handlers import exception_handler, sigint_handler
    from crisprme2.motifraptor_errors import Crisprme2FastaError
    import signal, os

    # Register the SIGINT handler at startup.
    signal.signal(signal.SIGINT, lambda sig, frame: sigint_handler())

    # Report a fatal I/O error from anywhere in the pipeline.
    try:
        with open(fastafile) as fh:
            data = fh.read()
    except IOError as e:
        exception_handler(
            Crisprme2FastaError,
            f"Failed reading Fasta file: {fastafile}",
            os.EX_IOERR,
            debug,
            e,
        )
"""

from typing import NoReturn, Optional
from colorama import init, Fore

import sys
import os


def sigint_handler() -> None:
    """Handle a SIGINT signal by printing a notice and exiting gracefully.

    Intended to be registered as the ``SIGINT`` signal handler at programme
    startup (e.g. via ``signal.signal(signal.SIGINT, lambda sig, frame:
    sigint_handler())``).  When a keyboard interrupt (``Ctrl+C``) is received,
    the function writes a brief notice to *stderr* and exits with
    ``os.EX_OSERR`` (exit code 71 on POSIX systems), which signals to the
    calling shell that the process was terminated by an OS-level event rather
    than a clean exit or a data error.

    Returns
    -------
    None
        This function never returns; it always terminates the process via
        :func:`sys.exit`.

    Raises
    ------
    SystemExit
        Always raised with ``os.EX_OSERR`` as the exit code.
    """
    # print message when SIGINT is caught to exit gracefully from the execution
    sys.stderr.write("\nCaught SIGINT. Exit CRISPRme2")
    sys.exit(os.EX_OSERR)  # mark as os error code


def exception_handler(
    exception_type: type,
    exception: str,
    code: int,
    debug: bool,
    e: Optional[Exception] = None,
) -> NoReturn:
    """Report a fatal error and either raise an exception or exit the process.

    This is the single point of error escalation for the entire CRISPRme2
    pipeline.  Its behaviour is controlled by the *debug* flag:

    * **Debug mode** (``debug=True``) — raises *exception_type* with a
      ``"\\n\\n"``-prefixed message so that the exception text is visually
      separated from the traceback in the terminal.  When *e* is provided,
      exception chaining (``raise ... from e``) is used to preserve the
      originating exception and its traceback.

    * **Normal mode** (``debug=False``) — writes a red ``ERROR:`` line to
      *stderr* (using ``colorama`` for cross-platform ANSI colour) and exits
      the process with *code*.  No traceback is printed.

    ``colorama.init()`` is called unconditionally at the start of every
    invocation.  This is safe because ``init()`` is idempotent; calling it
    multiple times has no side effects.

    Parameters
    ----------

    exception_type : type
        The exception class to instantiate and raise in debug mode.  Must be
        a subclass of :class:`~crisprme2.crisprme2_error.CrisprmeError`
        (or any :class:`Exception` subclass) and accept a single string
        argument.
    exception : str
        Human-readable description of the error.  In debug mode this string
        is passed to the exception constructor (with a leading ``"\\n\\n"``
        for readability); in normal mode it is written to *stderr*.
    code : int
        POSIX exit code used when terminating in normal mode.  Callers
        conventionally use ``os.EX_IOERR`` (74) for file I/O failures and
        ``os.EX_DATAERR`` (65) for data-parsing or data-consistency errors.
        Other ``os.EX_*`` constants may be used as appropriate.
    debug : bool
        When ``True``, raise *exception_type* with the full traceback intact.
        When ``False``, print a coloured error message and exit with *code*.
    e : Optional[Exception]
        Optional originating exception.  When provided in debug mode, it is
        used as the cause of the raised exception via ``raise ... from e``,
        preserving the full exception chain for diagnosis.  Ignored in normal
        mode.

    Returns
    -------
    NoReturn
        This function never returns; it either raises an exception or calls
        :func:`sys.exit`.

    Raises
    ------
    exception_type
        Raised in debug mode with *exception* as the message.  If *e* is not
        ``None``, the exception is chained from *e*.

    SystemExit
        Raised in normal mode with *code* as the exit status.

    Examples
    --------
    Report an I/O failure and exit (normal mode)::

        exception_handler(
            Crisprme2FastaError,
            f"Failed reading Fasta file: {fastafile}",
            os.EX_IOERR,
            debug=False,
        )

    Raise a chained exception with full traceback (debug mode)::

        try:
            df = pd.read_csv(fastafile, sep="\\t")
        except Exception as e:
            exception_handler(
                Crisprme2FastaError,
                f"Failed retrieving sequences from {mapfile}",
                os.EX_IOERR,
                debug=True,
                e=e,
            )
    """
    init()  # initialize colorama render
    if debug:  # debug mode -> always trace back the full error stack
        if e:  # inherits from previous error
            raise exception_type(f"\n\n{exception}") from e
        raise exception_type(f"\n\n{exception}")  # divide exception message from stack
    # gracefully trigger error and exit execution
    sys.stderr.write(f"{Fore.RED}\n\nERROR: {exception}\n{Fore.RESET}")
    sys.exit(code)  # exit execution returning appropriate error code
