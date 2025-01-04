#!/bin/bash
ROUNDS=5
ROOT_DIR=$(pwd)
LOG_ROOT="${ROOT_DIR}/result"
APP_ROOT="${ROOT_DIR}/system/ligra/apps"
DATA_ROOT="${ROOT_DIR}/dataset"

cd ${APP_ROOT} || exit
make clean && make -j
cd ${ROOT_DIR}

dataset=("MPI" "TW" "KR") #   
BFS_STARTER=("18767213" "3841886" "49100819") # "6557208"
SSSP_STARTER=("18767213" "3841886" "49100819") # "15620566"

for index in "${!dataset[@]}"
do 
    DATA=${dataset[index]}
    DATA_DIR=${DATA_ROOT}/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "=========="
    echo "dataset = ${DATA}"
    DATA_LOG_DIR=${LOG_ROOT}/ligra/${DATA}
    if [ ! -d "$DATA_LOG_DIR"  ]; then
        mkdir -p ${DATA_LOG_DIR}
    fi
    for APP in BFS BellmanFord Components PageRankDelta; do 
        echo "-----"
        echo "APP: ${APP}"
        APP_LOG=${DATA_LOG_DIR}/${APP}
        if [ ! -d "$APP_LOG"  ]; then
            mkdir -p ${APP_LOG}
        fi
        START_VERTEX_OPT=""
        ITER_OPT=""
        if [ "${APP}" == "BFS" ]; then
            START_VERTEX_OPT=" -r ${BFS_STARTER[index]} "
        fi
        if [ "${APP}" == "BellmanFord" ]; then
            START_VERTEX_OPT=" -r ${SSSP_STARTER[index]} "
        fi
        if [ "${APP}" == "PageRankDelta" ]; then
            ITER_OPT=" -iter 10 "
        fi 
        #   
        for ALGO in OG FO SO FBC HC DBG CO RBT GO; do # 
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
                echo "not exist reordering file : ${REORDER_FILE}"
                continue
            fi
            ${APP_ROOT}/${APP}  -rounds ${ROUNDS} ${START_VERTEX_OPT} \
            -o ${REORDER_FILE}  ${DATA_DIR}/${DATA}-${ALGO}.csr > ${LOG_FILE}
        done
    done
done

ext_data=("DL" "FR" "RM") #  
BFS_STARTER=("44365" "2681889" "55875322") # 
SSSP_STARTER=("44365" "2681889" "55875322") # 
DATA_ROOT="${ROOT_DIR}/ext_data"

for index in "${!ext_data[@]}"
do 
    DATA=${ext_data[index]}
    DATA_DIR=${DATA_ROOT}/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "=========="
    echo "dataset = ${DATA}"
    DATA_LOG_DIR=${LOG_ROOT}/ligra/${DATA}
    if [ ! -d "$DATA_LOG_DIR"  ]; then
        mkdir -p ${DATA_LOG_DIR}
    fi
    
    for APP in BFS BellmanFord Components PageRankDelta; do 
        echo "-----"
        echo "APP: ${APP}"
        APP_LOG=${DATA_LOG_DIR}/${APP}
        if [ ! -d "$APP_LOG"  ]; then
            mkdir -p ${APP_LOG}
        fi
        START_VERTEX_OPT=""
        ITER_OPT=""
        if [ "${APP}" == "BFS" ]; then
            START_VERTEX_OPT=" -r ${BFS_STARTER[index]} "
        fi
        if [ "${APP}" == "BellmanFord" ]; then
            START_VERTEX_OPT=" -r ${SSSP_STARTER[index]} "
        fi
        if [ "${APP}" == "PageRankDelta" ]; then
            ITER_OPT=" -iter 10 "
        fi 
        #   
        for ALGO in OG FO SO FBC HC DBG CO RBT GO; do # 
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
                echo "not exist reordering file : ${REORDER_FILE}"
                continue
            fi
            ${APP_ROOT}/${APP}  -rounds ${ROUNDS} ${START_VERTEX_OPT} \
            -o ${REORDER_FILE}  ${DATA_DIR}/${DATA}-${ALGO}.csr > ${LOG_FILE}
        done
    done
done
