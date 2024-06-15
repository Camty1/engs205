#!/bin/bash
cases=(4 8 16 32 64 128 256 512 1024)

for case in ${cases[@]}; do
    echo "N = $case"
    ./hw4 $case
done
