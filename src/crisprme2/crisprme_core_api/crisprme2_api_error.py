"""
crisprme2_api_error.py
---------------------
Custom exception hierarchy for the crisprme2 core API.

All public exceptions derive from Crisprme2Error so callers can catch the
entire family with a single ``except Crisprme2Error`` clause while still
being able to discriminate between subsystems when needed.

Hierarchy
~~~~~~~~~
::

    Crisprme2Error
    ├── Crisprme2BatcherError       - TargetBatcher failures
    ├── Crisprme2AnnotationError    - FeatureRegistry / annotation failures
    ├── Crisprme2AlignmentBatchError     - Alignment scoring / validation failures
    └── Crisprme2PipelineError      - Pipeline construction / runtime failures
            ├── Crisprme2PipelineConfigError   - bad arguments at build time
            ├── Crisprme2PipelineSubmitError   - failure while submitting a batch
            └── Crisprme2PipelineLifecycleError - close / context-manager misuse
"""

from ..crisprme2_error import Crisprme2Error


class Crisprme2BatcherError(Crisprme2Error):
    """Raised when the Rust TargetBatcher encounters an unrecoverable error"""

    def __init__(self, value: str) -> None:
        super().__init__(value)

    def __str__(self):
        return super().__str__()


class Crisprme2AnnotationError(Crisprme2Error):
    """Raised when feature registration or BED annotation fails"""

    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2AlignmentBatchError(Crisprme2Error):
    """Raised when an alignment attribute receives an invalid value"""

    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class Crisprme2PipelineError(Crisprme2Error):
    """
    Base class for all pipeline-related errors.

    Covers the full lifecycle of a :class:`~crisprme2.crisprme_core_api.Pipeline`:
    construction, batch submission, and shutdown.
    """

    def __init__(self, value: str) -> None:
        super().__init__(value)

    def __str__(self):
        return super().__str__()


class Crisprme2PipelineConfigError(Crisprme2PipelineError):
    """
    Raised when the pipeline cannot be constructed due to invalid arguments.

    Typical causes:
    - ``chunks`` is not a positive integer.
    - ``transforms`` is empty, contains non-callable items, or items that
      lack the required ``__call__`` signature.
    - The ``Thresholds`` object is missing or malformed.
    - The native ``pipeline()`` Rust function raises during initialisation.
    """

    def __init__(self, value: str) -> None:
        super().__init__(value)

    def __str__(self) -> str:
        return super().__str__()


class Crisprme2PipelineSubmitError(Crisprme2PipelineError):
    """
    Raised when submitting a :class:`~crisprme2.crisprme_core_api.TargetBatcher`
    batch to the pipeline fails.

    Typical causes:
    - The pipeline has already been closed.
    - The ``submit()`` call raises (e.g. channel disconnected, OOM).
    - The supplied ``batcher`` argument is not a ``TargetBatcher`` instance.
    """

    def __init__(self, value: str) -> None:
        super().__init__(value)

    def __str__(self) -> str:
        return super().__str__()


class Crisprme2PipelineLifecycleError(Crisprme2PipelineError):
    """
    Raised when the pipeline is used outside its valid lifecycle window.

    Typical causes:
    - Calling :meth:`~Pipeline.submit` after the context manager has exited.
    - Re-entering a pipeline context manager that has already been closed.
    - The Rust ``close()`` call raises during ``__exit__``.
    """

    def __init__(self, value: str) -> None:
        super().__init__(value)

    def __str__(self) -> str:
        return super().__str__()
