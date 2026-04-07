import crisprme2._crisprme2_native as nat
import numpy as np
import time

class Scorer:
    def __init__(self, score: int, value: float):
        self.score = score
        self.value = value

    def __call__(self, batch: nat.PyAlignmentBatch):
        scores = np.asarray(batch.score(self.score))
        scores[:] = self.value


class Printer:
    def __call__(self, batch: nat.PyAlignmentBatch):
        print("--- batch ---------------------------\n")

        self.print_simple_debug("seq_id", batch.seq_id())
        self.print_simple_debug("offset", batch.offset())

        self.print_str_debug("rguide", 32, batch.rguide())
        self.print_str_debug("rseq",   32, batch.rseq())

        self.print_simple_debug("score0", batch.score(0))
        self.print_simple_debug("score1", batch.score(1))
        self.print_simple_debug("score2", batch.score(2))
        self.print_simple_debug("score3", batch.score(3))

    def print_str_debug(self, name: str, len: int, view):
        array = np.asarray(view)
        print(f"+ {name} {array.shape}\n")
        for s in array.view(f'S{len}'):
            print(f"  {str(s)}")
        print()

    def print_simple_debug(self, name: str, view):
        array = np.asarray(view)
        print(f"+ {name} {array.shape}\n")
        print(" ", array)
        print()


nat.init_tracing()
pipeline = nat.dataset_pipeline(
    folder = "examples/data",
    batch_size = 10_000,
    guide = nat.Guide("GATTACA"),
    thresholds = nat.Thresholds(3, 3, 3),
    sequence_len = 32,
    chunks = 10_000, # 6 GB
    transforms = [
        Scorer(0, 2),
        Scorer(1, 7),
        #Printer(),
    ],
)