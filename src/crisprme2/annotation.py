""" """

from .crisprme_core_api import FeatureRegistry
from .bedfile import AnnotationBed

from typing import List, Tuple


# NOTE: targets_batch will be a list of views for each alignment to annotate
#           the views will point to contig, start/stop position of the target
def annotate(
    annotation: AnnotationBed,
    registry: FeatureRegistry,
    targets_batch: List[Tuple[str, int, int]],
) -> List[bytes]:
    return registry.annotate_batch(annotation, targets_batch)
