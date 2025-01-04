#!/bin/bash
ROUNDS=5
ROOT_DIR=$(pwd)
LOG_ROOT="${ROOT_DIR}/result"
GPOP_ROOT="${ROOT_DIR}/system/GPOP"
DATA_ROOT="${ROOT_DIR}/dataset"

dataset=("MPI" "TW" "KR") #   "TW" "KR"
BFS_STARTER=("18767213" "3841886" "49100819") # "6557208"
SSSP_STARTER=("18767213" "3841886" "49100819") # "15620566"

for index in "${!dataset[@]}"
do 
    DATA=${dataset[index]}
    DATA_DIR=${DATA_ROOT}/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "=========="
    echo "dataset = ${DATA}"
    DATA_LOG_DIR=${LOG_ROOT}/GPOP/${DATA}
    if [ ! -d "$DATA_LOG_DIR"  ]; then
        mkdir -p ${DATA_LOG_DIR}
    fi
    
    for APP in bfs sssp cc pr; do 
        APP_ROOT=${GPOP_ROOT}/${APP}
        cd ${APP_ROOT} || exit
        make clean && make -j
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
        for ALGO in OG FO SO FBC HC DBG CO RBT GO; do # 
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
                echo "not exist reordering file: ${REORDER_FILE}"
                continue
            fi
            ${APP_ROOT}/${APP} ${DATA_DIR}/${DATA}-${ALGO}.csr -rounds ${ROUNDS} \
            -map ${REORDER_FILE} ${START_VERTEX_OPT} ${ITER_OPT} > ${LOG_FILE}
        done
    done
done

ext_data=("DL" "FR" "RM") #   "FR" "RM"
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
    DATA_LOG_DIR=${LOG_ROOT}/GPOP/${DATA}
    if [ ! -d "$DATA_LOG_DIR"  ]; then
        mkdir -p ${DATA_LOG_DIR}
    fi
    for APP in bfs sssp cc pr; do 
        APP_ROOT=${GPOP_ROOT}/${APP}
        cd ${APP_ROOT} || exit
        make clean && make -j
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
        for ALGO in OG FO SO FBC HC DBG CO RBT GO; do # OG FO SO FBC HC DBG CO RBT GO
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
                echo "not exist reordering file: ${REORDER_FILE}"
                continue
            fi
            ${APP_ROOT}/${APP} ${DATA_DIR}/${DATA}-${ALGO}.csr -rounds ${ROUNDS} \
            -map ${REORDER_FILE} ${START_VERTEX_OPT} ${ITER_OPT} > ${LOG_FILE}
        done
    done
done

cd ${ROOT_DIR}
