#!/bin/bash
dimensions=(8 16 32 64 128 256)

for d in ${dimensions[@]}; do
    echo "N = $d"
    echo "Pseudospectral:"
    ./hw5_pseudo $d
    echo "Collocation:"
    ./hw5_coll $d
    echo "Antialiasing:"
    ./hw5_anti $d
done
