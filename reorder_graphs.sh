#!/bin/bash
PART_SIZE=512
FEAT_SIZE=10
ALGO=("FO" "FON" "SO" "FBC" "HC" "DBG" "CO") # 
make clean
make
for DATA in MPI TW KR;do 
    DATA_ROOT="dataset"
    DATA_DIR=${DATA_ROOT}/${DATA}
    FEAT_FILE=${DATA_DIR}/feat/bfs-15.feat
    REORDER_DIR=${DATA_DIR}/reorder
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
        if [ -f "$CSR_FILE" ];then
            rm ${CSR_FILE}
        fi
        ./frontorder -d ${DATA_DIR}/${DATA}.csr -a ${index} \
        -s ${PART_SIZE} -o ${CSR_FILE} -r $REORDER_FILE \
        -i ${FEAT_FILE} -k ${KV} 
    done
done

for DATA in DL FR RM;do 
    DATA_ROOT="ext_data"
    DATA_DIR=${DATA_ROOT}/${DATA}
    FEAT_FILE=${DATA_DIR}/feat/bfs-15.feat
    REORDER_DIR=${DATA_DIR}/reorder
    KV=4
    if [ "${DATA}" == "DL" ]; then
        KV=64
    fi
    if [ "${DATA}" == "FR" ]; then
        KV=32
    fi
    if [ "${DATA}" == "RM" ]; then
        KV=20
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
        if [ -f "$CSR_FILE" ];then
            rm ${CSR_FILE}
        fi
        ./frontorder -d ${DATA_DIR}/${DATA}.csr -a ${index} \
        -s ${PART_SIZE} -o ${CSR_FILE} -r $REORDER_FILE \
        -i ${FEAT_FILE} -k ${KV} 
    done
done
