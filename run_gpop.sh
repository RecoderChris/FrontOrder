#!/bin/bash
ROOT_DIR=$(pwd)
LOG_ROOT="${ROOT_DIR}/result"
GPOP_ROOT="${ROOT_DIR}/system/GPOP"
DATA_ROOT="${ROOT_DIR}/dataset"

ALGO=("FO" "SO" "FBC" "HC" "DBG" "CO")
dataset=("MPI" "TW" "KR") #  "DL"
BFS_STARTER=("16580161" "13243450" "18569402") # "6557208" 
SSSP_STARTER=("5655958" "13886266" "31631843") # "15620566"
ROUNDS=10

for index in "${!dataset[@]}"
do 
    DATA=${dataset[index]}
    DATA_DIR=${DATA_ROOT}/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "=========="
    echo "dataset = ${DATA}"
    if [ ! -d "$REORDER_DIR"  ];then
        echo "No reordering file!"
        continue
    fi
    DATA_LOG_DIR=${LOG_ROOT}/GPOP/${DATA}
    if [ ! -d "$DATA_LOG_DIR"  ]; then
        mkdir -p ${DATA_LOG_DIR}
    fi
    
    for APP in bfs sssp cc pr; do 
        APP_ROOT=${GPOP_ROOT}/${APP}
        cd ${APP_ROOT} || exit
        make clean
        make -j
        echo "-----"
        echo "APP: ${APP}"
        APP_LOG=${DATA_LOG_DIR}/${APP}
        if [ ! -d "$APP_LOG"  ]; then
            mkdir -p ${APP_LOG}
        fi
        START_VERTEX_OPT=""
        if [ "${APP}" == "bfs" ]; then
            START_VERTEX_OPT=" -s ${BFS_STARTER[index]} "
        fi
        if [ "${APP}" == "sssp" ]; then
            START_VERTEX_OPT=" -s ${SSSP_STARTER[index]} "
        fi
        if [ "${APP}" == "cc" ]; then
            NUMA_OPT="numactl -i all "
        fi
        if [ "${APP}" == "pr" ]; then
            ITER_OPT=" -iter 10 "
        fi 
        #   
        for ALGO in OG FO SO FBC HC DBG CO; do # 
            echo "------"
            echo "REORDER_ALGO: ${ALGO}"
            LOG_FILE=${APP_LOG}/${ALGO}.log
            if [ "${ALGO}" == "OG" ]; then
                ${APP_ROOT}/${APP} ${DATA_DIR}/${DATA}.csr \
                -rounds ${ROUNDS} ${START_VERTEX_OPT}  ${ITER_OPT} > ${LOG_FILE}
                continue
            fi
            REORDER_FILE=${REORDER_DIR}/${DATA}-${ALGO}.reorder
            if [ ! -f "$REORDER_FILE"  ];then
                echo ${REORDER_FILE}
                echo "not exist reordering file"
                continue
            fi
            ${APP_ROOT}/${APP} ${DATA_DIR}/${DATA}-${ALGO}.csr -rounds ${ROUNDS} \
            -map ${REORDER_FILE} ${START_VERTEX_OPT} ${ITER_OPT} > ${LOG_FILE}
        done
    done
done