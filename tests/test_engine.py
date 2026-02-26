
import crisprme2._crisprme2_native as native
import crisprme2

import pytest
import time

def setup_module():
    native.initialize_engine_logger()

def test_create_alignment_params():
    native.AlignmentParams(
        sequence_len=30,
        sequence_batch_size=1000,
        alignment_batch_size=100,
        thresholds=native.Thresholds(qgap=1, tgap=1, mism=1),
        mutation_max=100,
        guide=native.Guide("NNNNNN")
    )

def test_create_thresholds():
    native.Thresholds(qgap=1, tgap=1, mism=1)

def test_create_engine_1():
    native.HybridEngine(native.AlignmentParams(
        sequence_len=30,
        sequence_batch_size=1000,
        alignment_batch_size=100,
        thresholds=native.Thresholds(qgap=1, tgap=1, mism=1),
        mutation_max=100,
        guide=native.Guide("NNNNNN")
    ))

# Check if there is some global state problem at initialization of the engine
def test_create_engine_second_time():
    native.HybridEngine(native.AlignmentParams(
        sequence_len=30,
        sequence_batch_size=1000,
        alignment_batch_size=100,
        thresholds=native.Thresholds(qgap=1, tgap=1, mism=1),
        mutation_max=100,
        guide=native.Guide("NNNNNN")
    ))

def test_target_batcher():
    native.TargetBatcher(
        pam_seq="NNN",
        guide_seq="NNNNNNNNNNN",
        size=30,
        right=True,
        threads=4,
        batch_hits=40,
        max_unique=10,
        overlap_left=40
    )

def test_pipeline():
    
    batcher = native.TargetBatcher(
        pam_seq="NNN",
        guide_seq="NNNNNNNNNNN",
        size=30,
        right=True,
        threads=4,
        batch_hits=40,
        max_unique=10,
        overlap_left=40
    )

    engine = native.HybridEngine(native.AlignmentParams(
        sequence_len=30,
        sequence_batch_size=1000,
        alignment_batch_size=100,
        thresholds=native.Thresholds(qgap=1, tgap=1, mism=1),
        mutation_max=100,
        guide=native.Guide("NNNNNN")
    ))

    def check_alignment_result(result):
        assert isinstance(result, native.AlignmentBatchView)
        assert result.batcher_id() == batcher.id

        view = memoryview(result)

        # Check memory layout of the result
        assert view.format == 'QBxxxxxxxIBBxx'
        assert view.shape == (result.size(),)
        assert view.strides == (24,) # 24 bytes between alignments

    # First send and receive 
    engine.send(batcher)
    time.sleep(2)
    for result in engine.receive_blocking(batcher):
        check_alignment_result(result)

    # Second send and receive to check if the engine can be reused
    engine.send(batcher)
    for result in engine.receive_blocking(batcher):
        check_alignment_result(result)