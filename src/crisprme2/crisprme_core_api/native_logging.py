"""
native_logging.py
-----------------
Python wrapper around the Rust logging bridge exposed via PyO3.

The native core (``_crisprme2_native``) emits structured diagnostics through
the Rust ``tracing`` framework.  Calling :func:`init_native_logging` installs
a bridge that forwards every such event into the three CRISPRme2 Python
loggers, so Rust-side and Python-side diagnostics land in the *same* log
files:

=============  =====================================  =========================
tracing level  Python destination                     file(s)
=============  =====================================  =========================
``ERROR``      ``errorlog``                           ``errors.log``
``WARN``       ``basiclog`` + ``verboselog`` (INFO)   ``basic.log`` + verbose
``INFO``       ``basiclog`` + ``verboselog``          ``basic.log`` + verbose
``DEBUG``      ``verboselog``                         ``verbose.log``
=============  =====================================  =========================

The bridge is a *process-global* subscriber: it must be installed **once**,
early in the run, before any native subsystem (``TargetBatcher``,
``Pipeline``, ...) is constructed — otherwise the events emitted during their
start-up are not captured.

Typical usage
~~~~~~~~~~~~~
::

    from crisprme2.logger import CrisprmeLoggers
    from crisprme2.crisprme_core_api import init_native_logging

    loggers = CrisprmeLoggers(outdir)
    init_native_logging(loggers)          # Rust events now flow into the logs

    # ... build TargetBatcher / Pipeline / etc.; their native traces are
    #     mirrored into basic.log / verbose.log / errors.log

Notes
~~~~~
- **Idempotent**: the first call wins.  Subsequent calls are no-ops (a debug
  line is emitted) and do **not** re-point the bridge at a different
  ``CrisprmeLoggers`` bundle.
- **One-way**: there is no uninstall.  This mirrors the Rust side, where the
  global ``tracing`` subscriber is set exactly once via ``try_init``.
- ``WARN`` is remapped onto ``INFO`` (with a ``"WARNING: "`` textual marker)
  because the Python handlers use exact-level filters and expose no WARNING
  sink.  ``TRACE`` is dropped on the Rust side to keep hot-loop events off
  the interpreter.
"""

from __future__ import annotations

import os

from .crisprme2_api_error import Crisprme2LoggingError
from ..logger import CrisprmeLoggers

try:
    from .._crisprme2_native import init_logging as _rust_init_logging
except ImportError:
    # fallback for development/testing before the extension is built
    _rust_init_logging = None


# ==============================================================================
# module state
# ==============================================================================

#: The Rust subscriber is process-global and can only be set once; track that
#: here so repeated calls give clear feedback instead of silently doing nothing.
_installed: bool = False


# ==============================================================================
# internal helpers
# ==============================================================================


def _require_native(loggers: CrisprmeLoggers) -> None:
    """Raise if the native extension has not been compiled"""
    if _rust_init_logging is None:
        loggers.errorlog.log_raise_exception(
            "Rust init_logging function not exposed to Python. Ensure the "
            "native extension (_crisprme2_native) is compiled and installed.",
            os.EX_CANTCREAT,
            Crisprme2LoggingError,
        )


# ==============================================================================
# public wrapper
# ==============================================================================


def init_native_logging(loggers: CrisprmeLoggers) -> None:
    """
    Install the Rust -> Python logging bridge.

    Wires the native ``tracing`` subscriber so that events emitted anywhere
    in the Rust core are forwarded into *loggers* and written to the usual
    CRISPRme2 log files.  Call this **once**, as early as possible, before
    any native object is created.

    Parameters
    ----------
    loggers : CrisprmeLoggers
        The logger bundle whose ``basiclog`` / ``verboselog`` / ``errorlog``
        receive the forwarded Rust events.

    Raises
    ------
    Crisprme2LoggingError
        If *loggers* is not a :class:`~crisprme2.logger.CrisprmeLoggers`,
        if the native extension is unavailable, or if the Rust installer
        raises.

    Notes
    -----
    Idempotent: only the first call installs the bridge; later calls emit a
    debug line and return without changing the active loggers.
    """
    global _installed

    # Validate the bundle first: if it is malformed we cannot route the error
    # through it, so raise the typed exception directly.
    if not isinstance(loggers, CrisprmeLoggers):
        raise Crisprme2LoggingError(
            f"'loggers' must be a CrisprmeLoggers instance, "
            f"got {type(loggers).__name__!r}"
        )

    _require_native(loggers)  # ensure native rust api is installed

    if _installed:
        loggers.verboselog.debug(
            "Native logging bridge already installed; call ignored"
        )
        return

    loggers.verboselog.debug("Installing native (Rust) -> Python logging bridge")
    try:
        _rust_init_logging(loggers)  # type: ignore[misc]
    except Exception as e:
        loggers.errorlog.log_raise_exception(
            f"Failed to install native logging bridge: {e}",
            os.EX_UNAVAILABLE,
            Crisprme2LoggingError,
        )

    _installed = True
    loggers.basiclog.info(
        "Native logging bridge installed - Rust events routed to CRISPRme2 loggers"
    )
