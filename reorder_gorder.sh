#!/bin/bash

cd baseline/gorder
make clean
make -j
cd ../.. || exit

DATA_ROOT=dataset
EXT_DATA_ROOT=ext_data
GORDER_ROOT=baseline/gorder

for DATA in KR TF TW; do 
    DATA_DIR=${DATA_ROOT}/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "=========="
    echo "dataset = ${DATA}"
    if [ ! -d "$REORDER_DIR"  ];then
        mkdir -p ${REORDER_DIR}
    fi
    REORDER_FILE=${REORDER_DIR}/${DATA}-GO.reorder
    if [ -f "$REORDER_FILE"  ];then
        rm ${REORDER_FILE}
    fi
    ${GORDER_ROOT}/Gorder ${DATA_DIR}/${DATA}.el -o ${REORDER_FILE}
    CSR_FILE=${DATA_DIR}/${DATA}-GO.csr
    make clean
    make
    ./frontorder -d ${DATA_DIR}/${DATA}.csr -a 6 -o ${CSR_FILE} -i ${REORDER_FILE}
done

for DATA in DL FR RM; do 
    DATA_DIR=${EXT_DATA_ROOT}/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "=========="
    echo "dataset = ${DATA}"
    if [ ! -d "$REORDER_DIR"  ];then
        mkdir -p ${REORDER_DIR}
        echo "Reorder Folder created."
    fi
    REORDER_FILE=${REORDER_DIR}/${DATA}-GO.reorder
    if [ -f "$REORDER_FILE"  ];then
        rm ${REORDER_FILE}
    fi
    ${GORDER_ROOT}/Gorder ${DATA_DIR}/${DATA}.el -o ${REORDER_FILE} 
    CSR_FILE=${DATA_DIR}/${DATA}-GO.csr
    make clean
    make
    ./frontorder -d ${DATA_DIR}/${DATA}.csr -a 6 -o ${CSR_FILE} -i ${REORDER_FILE}
done
