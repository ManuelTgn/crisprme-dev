#!/bin/bash

# Positive strand
./target/release/crisprme           \
  --sequence-batch-size 1000000     \
  --alignment-batch-size 10000000   \
  mine                              \
  --qgap 1 --tgap 1 --mism 4        \
  sequences_ref_chr1_NGG_positive   \
  22 CTAACAGTTGCTTTTATCAC \
  mine_chr1_positive.bin

# Negative strand
#./target/release/crisprme           \
#  --sequence-batch-size 1000000     \
#  --alignment-batch-size 10000000   \
#  mine                              \
#  --qgap 1 --tgap 1 --mism 4        \
#  sequences_ref_chr1_NGG_negative   \
#  22 CTAACAGTTGCTTTTATCAC \
#  mine_chr1_negative.bin
