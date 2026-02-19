""" """

from .bedfile import AnnotationBed, AnnotationGff

from typing import Union


def annotate_target(
    annotation_file: Union[AnnotationBed, AnnotationGff],
    contig: str,
    start: int,
    stop: int,
) -> str:
    annotation = annotation_file.fetch_features(contig, start, stop)
    if annotation is None:
        return "NA"
    return annotation
