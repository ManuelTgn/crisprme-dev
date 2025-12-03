#!/usr/bin/env bash

input_folder=$1
for f in $1/*; do
	echo "PROCESSING $f ..."
	./target/release/crisprme preprocess-list $f 30
	echo "DONE"
done
