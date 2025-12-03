
import os
import subprocess
import time

from multiprocessing import Pool, Queue, Process, cpu_count, current_process
from queue import Empty

# Worker for a single GPU
def worker(queue: Queue, gpu_id: int):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    while True:
        try:
            chrom = queue.get_nowait()
            command = f"./target/release/crisprme --sequence-batch-size 1000000 --alignment-batch-size 60000000 mine --qgap 1 --tgap 1 --mism 4 {chrom} 20 ATTGAGATAGTGGNGG mined.bin"
            print(f"command: {command}")
        except Empty:
            break

        start = time.time()
        subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        elapsed = time.time() - start

        print(f"gpu:{gpu_id}: complete for {chrom} in {elapsed} s")


CHROMS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY"
]

def main():

    queue = Queue()
    for chrom in CHROMS:
        queue.put(chrom)

    start = time.time()

    # spawn one worker per GPU
    workers = []
    for gpu_id in range(2):
        p = Process(target=worker, args=(queue, gpu_id))
        workers.append(p)
        p.start()

    # wait for all of them
    for p in workers:
        p.join()

    elapsed = time.time() - start
    print(f"benchmark complete in {elapsed} s")


if __name__ == "__main__":
    main()
