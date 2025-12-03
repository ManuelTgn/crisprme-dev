
import os
import time
import subprocess

from multiprocessing import Pool, Queue, Process, cpu_count, current_process
from queue import Empty

INPUT = [
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

SBS = 1000000
ABS = 1000000

def run(queue: Queue, gpu_id: int, guide, slen, ggap, sgap, mism):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    while True:
        try:
            chrom = queue.get_nowait()
            alignments = f"runs/alignments_{chrom}_{guide}_ggap{ggap}_sgap{sgap}_mism{mism}.bin"
            results = f"runs/results_{chrom}_{guide}_ggap{ggap}_sgap{sgap}_mism{mism}.csv"

            # Mine
            start = time.perf_counter()
            try:
                cmd = f"./target/release/crisprme --sequence-batch-size {SBS} --alignment-batch-size {ABS} mine --qgap {ggap} --tgap {sgap} --mism {mism} {chrom} {slen} {guide} {alignments}"
                print(f"GPU:{gpu_id}: command: {cmd}")
                subprocess.run(cmd, shell=True, check=True)

            except subprocess.CalledProcessError as e:
                print(f"GPU:{gpu_id}: failed with exit code {e.returncode}")
                return 0

            elapsed = time.perf_counter() - start
            print(f"GPU:{gpu_id}: elapsed time for {chrom}: {elapsed:.2f} s")

        except Empty:
            break


def spawn(gpu_count, guide, slen, ggap, sgap, mism):
    
    queue = Queue()
    for chrom in INPUT:
        queue.put(chrom)

    start = time.perf_counter()

    # Spawn one worker for each GPU
    workers = []
    for gpu_id in range(gpu_count):
        p = Process(target=run, args=(queue, gpu_id, guide, slen, ggap, sgap, mism))
        workers.append(p)
        p.start()

    # Wait for all runs to finish
    for p in workers:
        p.join()
    
    total_elapsed = time.perf_counter() - start
    print(f"total_elapsed: {total_elapsed:.2f} s")

    # Generate csv with results
    for chrom in INPUT:
        alignments = f"runs/alignments_{chrom}_{guide}_ggap{ggap}_sgap{sgap}_mism{mism}.bin"
        results = f"runs/results_{chrom}_{guide}_ggap{ggap}_sgap{sgap}_mism{mism}.csv"

        try:
            # crisprme results <INPUT> <ALIGNMENTS> <SEQUENCE_LEN> <GUIDE>
            cmd = f"./target/release/crisprme results {chrom} {alignments} {slen} {guide} --skip-wildcards > {results}"
            print(f"command: {cmd}")
            subprocess.run(cmd, shell=True, check=True)

        except subprocess.CalledProcessError as e:
            print(f"failed with exit code {e.returncode}")
            return 0

if __name__ == "__main__":
    GPU_COUNT = 2

    spawn(GPU_COUNT, "CTAACAGTTGCTTTTATCAC", 22, 1, 1, 4)
