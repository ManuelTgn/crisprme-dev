#!/usr/bin/env bash

# All chromosomes we want to preprocess
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

input=$1
for chr in "${CHROMS[@]}"; do
	echo "Processing $chr.fa ..."
	./target/release/crisprme preprocess "$input/$chr.fa" 22 1
done
