#!/bin/bash

make clean
make -j

DATA_ROOT=/home/zhangxm/24MLSys-HisOrder/data
RABBIT_ROOT=/home/zhangxm/24MLSys-HisOrder/baseline/rabbit
SRC_DIR=${RABBIT_ROOT}/demo

#for DATASET in web-Google usa-road twitter 
for DATA in KR ; do # web-Google usa-roa LJ DL TF FR OR WB TW UK KR AR UR UNI
    DATA_DIR=${DATA_ROOT}/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "==================="
    echo "[DATASET: ] ${DATA}"
    REORDER_FILE=${REORDER_DIR}/${DATA}-rabbit.reorder
    if [ -f "$REORDER_FILE"  ];then
      rm ${REORDER_FILE}
    fi
    if [ ! -d "$REORDER_DIR"  ];then
        mkdir -p ${REORDER_DIR}
        echo "Reorder Folder created."
    fi
    ${SRC_DIR}/reorder ${DATA_DIR}/${DATA}.el >> ${REORDER_FILE}
done
