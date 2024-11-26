#!/bin/bash

PART_SIZE=512
FEAT_SIZE=15
DATA_ROOT="dataset"

ALGO=("FO" "SO" "FBC" "HC" "DBG" "CO")

for DATA in MPI TW KR;do 
    DATA_DIR=${DATA_ROOT}/${DATA}
    FEAT_FILE=${DATA_DIR}/feat/bfs-15.feat
    REORDER_DIR=${DATA_DIR}/reorder
    KV=4
    if [ "${DATA}" == "DL" ]; then
        KV=64
    fi
    if [ "${DATA}" == "MPI" ]; then
        KV=16
    fi
    if [ "${DATA}" == "TW" ]; then
        KV=32
    fi
    if [ "${DATA}" == "KR" ]; then
        KV=30
    fi
    echo "=========="
    echo "dataset = ${DATA}"
    if [ ! -d "$REORDER_DIR"  ];then
        mkdir -p ${REORDER_DIR}
    fi
    
    for index in "${!ALGO[@]}"
    do
        echo "algo = ${ALGO[index]}"
        REORDER_ALGO="${ALGO[index]}"
        REORDER_FILE=${REORDER_DIR}/${DATA}-${REORDER_ALGO}.reorder
        CSR_FILE=${DATA_DIR}/${DATA}-${REORDER_ALGO}.csr
        if [ -f "$REORDER_FILE" ];then
            rm ${REORDER_FILE}
        fi
        ./frontorder -d ${DATA_DIR}/${DATA}.csr -a ${index} \
        -s ${PART_SIZE} -o ${CSR_FILE} -r $REORDER_FILE \
        -i ${FEAT_FILE} -k ${KV} 
    done
done