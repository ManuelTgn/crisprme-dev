#!/bin/bash
./target/release/crisprme --sequence-batch-size 10000000 --alignment-batch-size 1000000 mine --qgap 1 --tgap 1 --mism 4 chr1 22 CTAACAGTTGCTTTTATCAC mined.bin
