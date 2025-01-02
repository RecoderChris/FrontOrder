#!/bin/bash
ROOT_DIR=$(pwd)
LOG_ROOT="${ROOT_DIR}/result"
APP_ROOT="${ROOT_DIR}/system/ligra/apps"
DATA_ROOT="${ROOT_DIR}/dataset"

ALGO=("FO" "SO" "FBC" "HC" "DBG" "CO")
dataset=("MPI" "TW" "KR") #   "DL"
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
    DATA_LOG_DIR=${LOG_ROOT}/ligra/${DATA}
    if [ ! -d "$DATA_LOG_DIR"  ]; then
        mkdir -p ${DATA_LOG_DIR}
    fi
    
    for APP in BFS BellmanFord Components PageRankDelta; do 
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
            START_VERTEX_OPT=" -r ${BFS_STARTER[index]} "
        fi
        if [ "${APP}" == "sssp" ]; then
            START_VERTEX_OPT=" -r ${SSSP_STARTER[index]} "
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
                ${APP_ROOT}/${APP}  -rounds ${ROUNDS} ${START_VERTEX_OPT} \
                ${ITER_OPT} ${DATA_DIR}/${DATA}.csr > ${LOG_FILE}
                continue
            fi
            REORDER_FILE=${REORDER_DIR}/${DATA}-${ALGO}.reorder
            if [ ! -f "$REORDER_FILE"  ];then
                echo ${REORDER_FILE}
                echo "not exist reordering file"
                continue
            fi
            ${APP_ROOT}/${APP}  -rounds ${ROUNDS} ${START_VERTEX_OPT} \
            -o ${REORDER_FILE}  ${DATA_DIR}/${DATA}-${ALGO}.csr > ${LOG_FILE}
        done
    done
done