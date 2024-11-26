#!/bin/bash
NUMA_OPT="numactl -i all "

GPOP_ROOT=/home/zhangxm/24MLSys-HisOrder/benchmark-systems/ligra-master/apps
DATA_ROOT=/home/zhangxm/24MLSys-HisOrder/data

# datasets=("DL" "TF" "FR" "TW" "KR" "RM27")
# starters=("44365" "18767213" "2681889" "3841886" "49100819" "55875322") # 
datasets=("R22-4") #   "R22-8" "R22-16" "R22-32" "R22-64" "R22-128"
# datasets=("R22-64") #  
starter=("0" "0"  "0" "0" "10" "0") #   
# starter=("10") #   

# 构建字符串数组
for ((i=0; i<${#datasets[@]}; i++)); do
    DATA=${datasets[i]}
    START=${starter[i]}
    DATA_DIR=${DATA_ROOT}/${DATA}
    # GPOP_LOG_DIR=${DATA_DIR}/GPOP
    GPOP_LOG_DIR=/home/zhangxm/24MLSys-HisOrder/result/ligra_rmat_scale/${DATA}
    REORDER_DIR=${DATA_DIR}/reorder
    echo "==================="
    echo "DATASET: ${DATA}"
    if [ ! -d "$GPOP_LOG_DIR"  ];then
        mkdir -p "${GPOP_LOG_DIR}"
        echo "ligra log dir created." 
    fi
    if [ ! -d "$REORDER_DIR"  ];then
        echo "No reordering file!"
        continue
    fi

    for APP in PageRank; do # BFS Components PageRank Components PageRankDelta BellmanFord
        cd ${GPOP_ROOT} || exit
        echo "-------------------"
        echo "APP: ${APP}"
        APP_LOG=${GPOP_LOG_DIR}/${APP}
        if [ ! -d "$APP_LOG"  ]; then
            mkdir -p ${APP_LOG}
            echo "${APP} Folder created."
        fi
        for ALGO in Orig HOSC; do # Orig HO-SC Rand gorder rabbit HubCluster HubSort DBG Sort Corder  Hisorder Hisorder_noblc HO-embed HN-embed HOSC-embed
            echo "-------------------"
            echo "REORDER_ALGO: ${ALGO}"
            LOG_FILE=${APP_LOG}/${ALGO}.log
            if [ -f "$LOG_FILE"  ];then
                rm ${LOG_FILE}
            fi
            if [ "${ALGO}" == "Orig" ]; then
                ${NUMA_OPT} ${GPOP_ROOT}/${APP} -r "${START}" -rounds 5  ${DATA_DIR}/${DATA}.csr >> "${LOG_FILE}"
                continue
            fi
            REORDER_FILE=${REORDER_DIR}/${DATA}-${ALGO}.reorder
            if [ ! -f "$REORDER_FILE"  ];then
                echo "No reordering file!"
                continue
            fi
            ${NUMA_OPT} ${GPOP_ROOT}/${APP} -r "${START}" -rounds 5 -o ${REORDER_FILE}  ${DATA_DIR}/${DATA}-${ALGO}.csr >> "${LOG_FILE}" 
        done
    done
done
