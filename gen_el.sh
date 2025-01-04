#!/bin/bash

# This script transforms CSR format(*.csr) to edgelist format
# The edgelist format is required by rabbit and gorder reordering methods
cd tools/gen_el || exit
make clean && make

ROOT=../../dataset
for DATA in MPI KR TW ; do # 
    DATA_DIR=${ROOT}/${DATA}
    echo "========="
    echo "Dataset : ${DATA}"
    # generate edgelist
    ./gen_el ${DATA_DIR}/${DATA}.csr ${DATA_DIR}/${DATA}.el
done

EXT_ROOT=../../ext_data
for DATA in DL FR RM; do # DL FR RM
    DATA_DIR=${EXT_ROOT}/${DATA}
    echo "========="
    echo "Dataset : ${DATA}"
    # generate edgelist
    ./gen_el ${DATA_DIR}/${DATA}.csr ${DATA_DIR}/${DATA}.el
done

cd ../..||exit
