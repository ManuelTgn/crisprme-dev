import subprocess
import time
import csv

from pathlib import Path

global_total_elapsed = 0

def benchmark(chrom, runs, output_file):
    results = []

    total_elapsed = 0
    for i in range(runs):
        start = time.perf_counter()
        try:
            #command = f"CUDA_VISIBLE_DEVICES=1 ./target/release/crisprme --sequence-batch-size 10000000 --alignment-batch-size 1000000 mine --qgap 1 --tgap 1 --mism 4 {chrom} 22 CTAACAGTTGCTTTTATCAC mined.bin"
            command = f"CUDA_VISIBLE_DEVICES=1 ./target/release/crisprme --sequence-batch-size 1000000 --alignment-batch-size 10000000 mine --qgap 1 --tgap 1 --mism 2 {chrom} 20 ATTGAGATAGTGGNGG mined.bin"
            print(f"command: {command}")
            subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Run {i+1} failed with exit code {e.returncode}")
            elapsed = None
        else:
            elapsed = time.perf_counter() - start
            total_elapsed += elapsed


        results.append((i + 1, elapsed))
        print(f"Run {i+1}: {elapsed:.6f} seconds" if elapsed is not None else f"Run {i+1}: FAILED")


    # Save to file (plain text)
    output_path = Path(output_file)
    with output_path.open("w") as f:
        f.write("Run,ElapsedTimeSeconds\n")
        for run, elapsed in results:
            f.write(f"{run},{elapsed if elapsed is not None else 'FAILED'}\n")

    return total_elapsed / runs


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

if __name__ == "__main__":
    for chrom in INPUT:
        global_total_elapsed += benchmark(chrom, 1, f"benchmark/{chrom}.csv")
        print(f"total elapsed from beginning: {global_total_elapsed}s")
