#!/bin/bash

cd baseline/rabbit/demo
make clean
make -j
cd ../../.. || exit

DATA_ROOT=dataset
EXT_DATA_ROOT=ext_data
RBT_ROOT=baseline/rabbit/demo

for DATA in KR TF TW; do 
    DATA_DIR=${DATA_ROOT}/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "=========="
    echo "dataset = ${DATA}"
    if [ ! -d "$REORDER_DIR"  ];then
        mkdir -p ${REORDER_DIR}
    fi
    REORDER_FILE=${REORDER_DIR}/${DATA}-RBT.reorder
    if [ -f "$REORDER_FILE"  ];then
        rm ${REORDER_FILE}
    fi
    
    ${RBT_ROOT}/reorder ${DATA_DIR}/${DATA}.el >> ${REORDER_FILE}

    CSR_FILE=${DATA_DIR}/${DATA}-RBT.csr
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
    fi
    REORDER_FILE=${REORDER_DIR}/${DATA}-RBT.reorder
    if [ -f "$REORDER_FILE"  ];then
        rm ${REORDER_FILE}
    fi
    
    ${RBT_ROOT}/reorder ${DATA_DIR}/${DATA}.el >> ${REORDER_FILE}
    CSR_FILE=${DATA_DIR}/${DATA}-RBT.csr
    make clean
    make
    ./frontorder -d ${DATA_DIR}/${DATA}.csr -a 6 -o ${CSR_FILE} -i ${REORDER_FILE}
done
