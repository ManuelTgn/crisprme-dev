
import crisprme2._crisprme2_native as native
import crisprme2

import numpy as np
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
    
    batcher_A = native.TargetBatcher(
        pam_seq="NNN",
        guide_seq="NNNNNNNNNNN",
        size=30,
        right=True,
        threads=4,
        batch_hits=40,
        max_unique=10,
        overlap_left=40
    )

    batcher_B = native.TargetBatcher(
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

    def check_alignment_result(result, batcher):
        assert isinstance(result, native.AlignmentBatchView)
        assert result.batcher_id() == batcher.id

        view = memoryview(result)
        print(view.shape, view.format, view.strides, view.nbytes)

        # Check memory layout of the result
        assert view.format == 'QBxxxxxxxIBBxx'
        assert view.shape == (result.size(),)
        assert view.strides == (24,) # 24 bytes between alignments

        dt = np.dtype([
            ('cigarx_storage', np.uint64, 1),
            ('cigarx_bits',    np.uint8,  1),
            ('_pad_0',         np.uint8,  7),
            ('id',             np.uint32, 1),
            ('offset',         np.uint8,  1),
            ('strand',         np.uint8,  1),
            ('_pad_1',         np.uint8,  2)
        ])

        np_arr = np.frombuffer(view, dtype=dt)
        assert np_arr.shape == (result.size(),)
        print(np_arr)

    # First send and receive 
    engine.send(batcher_B)
    for result in engine.receive_blocking(batcher_B):
        check_alignment_result(result, batcher_B)

    # Second send and receive to check if the engine can be reused
    engine.send(batcher_A)
    for result in engine.receive_blocking(batcher_A):
        check_alignment_result(result, batcher_A)