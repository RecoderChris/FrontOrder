#!/bin/bash

DATA_ROOT=/home/zhangxm/24MLSys-HisOrder/data

for DATA in KR; do #TW TF FR 
    DATA_DIR=${DATA_ROOT}/${DATA}
    for ALGO in Orig HisOrder Hisorder_noblc HubCluster; do
        DATA_EL=${DATA_DIR}/${DATA}-${ALGO}.el
        DATA_ADJ=${DATA_DIR}/${DATA}-${ALGO}.adj
        if [ "${ALGO}" == "Orig" ]; then
            DATA_EL=${DATA_DIR}/${DATA}.el
            DATA_ADJ=${DATA_DIR}/${DATA}.adj
        fi
        ./SNAPtoAdj ${DATA_EL} ${DATA_ADJ}
    done
done