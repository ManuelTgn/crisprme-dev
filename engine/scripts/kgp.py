
import os
import time
import subprocess

from multiprocessing import Pool, Queue, Process, cpu_count, current_process
from queue import Empty

INPUT = [
"kgp_chr10_copy0",
"kgp_chr10_copy1",
"kgp_chr11_copy0",
"kgp_chr11_copy1",
"kgp_chr12_copy0",
"kgp_chr12_copy1",
"kgp_chr13_copy0",
"kgp_chr13_copy1",
"kgp_chr14_copy0",
"kgp_chr14_copy1",
"kgp_chr15_copy0",
"kgp_chr15_copy1",
"kgp_chr16_copy0",
"kgp_chr16_copy1",
"kgp_chr17_copy0",
"kgp_chr17_copy1",
"kgp_chr18_copy0",
"kgp_chr18_copy1",
"kgp_chr19_copy0",
"kgp_chr19_copy1",
"kgp_chr1_copy0",
"kgp_chr1_copy1",
"kgp_chr20_copy0",
"kgp_chr20_copy1",
"kgp_chr21_copy0",
"kgp_chr21_copy1",
"kgp_chr22_copy0",
"kgp_chr22_copy1",
"kgp_chr2_copy0",
"kgp_chr2_copy1",
"kgp_chr3_copy0",
"kgp_chr3_copy1",
"kgp_chr4_copy0",
"kgp_chr4_copy1",
"kgp_chr5_copy0",
"kgp_chr5_copy1",
"kgp_chr6_copy0",
"kgp_chr6_copy1",
"kgp_chr7_copy0",
"kgp_chr7_copy1",
"kgp_chr8_copy0",
"kgp_chr8_copy1",
"kgp_chr9_copy0",
"kgp_chr9_copy1",
"kgp_chrX_copy0",
"kgp_chrX_copy1",
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

    spawn(GPU_COUNT, "CTAACAGTTGCTTTTATCAC", 30, 1, 1, 4)
