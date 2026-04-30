"""
callbacks.py
Structural interfaces and base implementations for alignment batch processing.

This module defines the AlignmentTransformer protocol, which serves as the
foundation for all batch-level modifications and inspections within the pipeline.
By utilizing Python's structural subtyping (Protocols), it ensures that
disparate tools—ranging from scorers to debug loggers—share a consistent
execution contract.

Typical usage

::

    from crisprme2.transforms import Scorer, Printer
    from crisprme2.callbacks import CallbackPipeline

    # Define a pipeline of operations
    pipeline = CallbackPipeline([
        Scorer(score_idx=0, value=1.0),
        Printer()
    ])

    # The pipeline itself satisfies the AlignmentTransformer protocol
    pipeline(batch)

Notes
~~~~~
- All transformers must implement the ``__call__`` method accepting a
  ``nat.PyAlignmentBatch``.
- Most transformers modify the batch **in-place** to maximize performance
  and minimize memory overhead in high-throughput bioinformatics tasks.
- The ``@runtime_checkable`` decorator is used to allow ``isinstance``
  checks against the protocol at runtime.
"""

from __future__ import annotations

from typing import Any, List, Optional, Protocol, runtime_checkable

import os

from .crisprme_core_api import AlignmentBatch
from .crisprme2_error import Crisprme2CallbackError
from .logger import CrisprmeLoggers

# ------------------------------------------------------------------------------
# public protocols
# ------------------------------------------------------------------------------


@runtime_checkable
class Transformer(Protocol):
    """
    Structural interface for alignment batch transformations.

    Any class implementing this protocol can be passed to the pipeline
    engine to inspect or modify alignment data. Implementations should
    handle the ``nat.PyAlignmentBatch`` object directly, often leveraging
    NumPy views for performance.
    """

    def __call__(self, batch: AlignmentBatch) -> None:
        """
        Execute the transformation on the provided batch.

        Parameters
        ----------
        batch : nat.PyAlignmentBatch
            The alignment batch object to be processed/modified.
        """
        ...


# ------------------------------------------------------------------------------
# composite implementations
# ------------------------------------------------------------------------------


class CallbackPipeline:
    """
    A composite container for sequencing multiple alignment transformers.

    This class aggregates multiple ``AlignmentTransformer`` instances and
    executes them sequentially on a single batch. It implements the
    transformer protocol itself, allowing for nested pipelines.

    Parameters
    ----------
    callbacks : List[AlignmentTransformer], optional
        A list of initialized transformer objects to execute in order.

    Attributes
    ----------
    callbacks : List[AlignmentTransformer]
        The internal list of registered callbacks.
    """

    def __init__(
        self, loggers: CrisprmeLoggers, callbacks: Optional[List[Transformer]] = None
    ) -> None:
        self._callbacks: List[Transformer] = callbacks or []
        self._loggers = loggers

    def add(self, callback: Transformer) -> None:
        """
        Append a new transformer to the end of the pipeline.

        Parameters
        ----------
        callback : AlignmentTransformer
            The transformer instance to register.
        """
        self._callbacks.append(callback)

    def __call__(self, batch: AlignmentBatch) -> Any:
        """
        Execute all registered callbacks on the batch.

        Callbacks are executed in the order they were added. If a callback
        fails, the error is logged and re-raised to stop the pipeline.

        Parameters
        ----------
        batch : nat.PyAlignmentBatch
            The alignment batch object to be processed.

        Raises
        ------
        Exception
            Any exception raised by an individual callback.
        """
        for callback in self._callbacks:
            try:
                callback(batch)
            except Exception as e:
                self._loggers.errorlog.log_raise_exception(
                    f"Callback execution failed ({callback.__call__.__name__}): {e}",
                    os.EX_PROTOCOL,
                    Crisprme2CallbackError,
                )

    # --------------------------------------------------------------------------
    # value-type protocol
    # --------------------------------------------------------------------------

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(callbacks={len(self._callbacks)})"
