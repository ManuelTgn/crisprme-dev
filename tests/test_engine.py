
import crisprme2._crisprme2_native as m
import time

def test_create_alignment_params():
    m.AlignmentParams(
        sequence_len=30,
        sequence_batch_size=1000,
        alignment_batch_size=100,
        thresholds=m.Thresholds(qgap=1, tgap=1, mism=1),
        mutation_max=100,
        guide=m.Guide("NNNNNN")
    )

def test_create_thresholds():
    m.Thresholds(qgap=1, tgap=1, mism=1)

def test_create_engine():
    m.HybridEngine(m.AlignmentParams(
        sequence_len=30,
        sequence_batch_size=1000,
        alignment_batch_size=100,
        thresholds=m.Thresholds(qgap=1, tgap=1, mism=1),
        mutation_max=100,
        guide=m.Guide("NNNNNN")
    ))

