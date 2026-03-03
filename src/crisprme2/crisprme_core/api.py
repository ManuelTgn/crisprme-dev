""" """

from ..bedfile import AnnotationBed
from ..logger import CrisprmeLoggers
from .crisprme_api_error import Crisprme2AnnotationError

from .._crisprme2_native import PyRegistry

from typing import List, Tuple

import os


class FeatureRegistry:

    def __init__(self, path: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        try:  # initialize rust-based feature registry
            self._registry = PyRegistry(path)
            self._num_features = self._registry.num_features()
        except Exception as e:
            loggers.errorlog.log_raise_exception(
                f"Feature registry initialization failed on {path}: {e}",
                os.EX_IOERR,
                Crisprme2AnnotationError,
            )

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
                _get_feature_id(self._registry, f, self._loggers) for f in features
            ]
            hits_per_target.append(features_ids)
        # rust parallel batch annotation
        return self._registry.annotate_batch(hits_per_target)

    @property
    def num_features(self) -> int:
        return self._num_features


def _get_feature_id(
    registry: PyRegistry, feature: str, loggers: CrisprmeLoggers
) -> int:
    # recover feature id in current registry
    fid = registry.get_feature_id(feature)
    if fid is None:
        loggers.errorlog.log_raise_exception(
            f"Feature '{feature}' not found in registry",
            os.EX_DATAERR,
            Crisprme2AnnotationError,
        )
    return fid
