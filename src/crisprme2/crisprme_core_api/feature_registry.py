""" """

from ..bedfile import AnnotationBed
from ..logger import CrisprmeLoggers
from .crisprme_api_error import Crisprme2AnnotationError, Crisprme2AlignmentError

try:  # import rust API modules
    from .._crisprme2_native import PyRegistry as RustFeatureRegistry
except ImportError:
    # fallback for development/testing
    RustFeatureRegistry = None

from typing import List, Tuple, Optional

import os


class FeatureRegistry:

    def __init__(self, path: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        if RustFeatureRegistry is None:
            self._loggers.errorlog.log_raise_exception(
                "Rust FeatureRegistry not exposed to python",
                os.EX_CANTCREAT,
                ValueError,
            )
        # load annotation features and register action on log
        self._loggers.basiclog.info(f"Loading features from {path}")
        try:  # initialize rust-based feature registry
            self._path = path
            self._registry = RustFeatureRegistry(path)
            self._num_features = self._registry.num_features()
        except Exception as e:
            loggers.errorlog.log_raise_exception(
                f"Feature registry initialization failed on {path}: {e}",
                os.EX_IOERR,
                Crisprme2AnnotationError,
            )
        self._loggers.verboselog.info(
            f"FeatureRegistry initialized with {self._num_features} unique features"
        )

    @property
    def path(self) -> str:
        return self._path

    @property
    def num_features(self) -> int:
        return self._num_features

    def _get_feature_id(self, feature_name: str) -> Optional[int]:
        return self._registry.get_feature_id(feature_name)

    def annotate_batch(
        self, annotation: AnnotationBed, targets: List[Tuple[str, int, int]]
    ) -> List[bytes]:
        if not targets:
            self._loggers.errorlog.log_raise_exception(
                "Empty targets batch", os.EX_DATAERR, Crisprme2AnnotationError
            )
        # collect feature id hits per target
        hits_per_target: List[List[int]] = []
        for contig, start, stop in targets:
            features = annotation.fetch_features(contig, start, stop)
            # no overlap between targets and features -> empty list
            if not features:
                hits_per_target.append([])  # must not be empty
                continue
            # convert feature names to ids using rust registry
            features_ids = [
                fid for f in features if (fid := self._get_feature_id(f)) is not None
            ]
            hits_per_target.append(features_ids)
        # rust parallel batch annotation
        return self._registry.annotate_batch(hits_per_target)
