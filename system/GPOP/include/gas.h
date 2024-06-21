/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * This code implements work optimized propagation blocking with
 * transposed bin graph to reduce cache misses in scatter
 */


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include <vector>
//#include <atomic>
#include "../include/partition.h"

#define SPARSE_MODE
// #undef SPARSE_MODE


using namespace std;


template<class graph>
void getFrontier(graph* G)
{
    intV* prefix = new intV [(G->partListPtr)+1]();
    prefix[0] = 0;
    for (intV i=0; i<G->partListPtr; i++)
    {
        prefix[i+1] = prefix[i] + G->TD[G->activeScatter[i]].frontierSize;
    } 
    #pragma omp parallel for
    for (intV i=0; i<G->partListPtr; i++)
    {
        for (intV j=0; j<G->TD[G->activeScatter[i]].frontierSize; j++)
            G->frontier[prefix[i]+j] = G->TD[G->activeScatter[i]].frontier[j];
    }
    delete[] prefix;
}


template<class graph>
void resetFrontier(graph* G)
{
    #pragma omp parallel for
    for (intV i=0; i<G->partListPtr; i++)
    {
        for (intV j=0; j<G->TD[G->activeScatter[i]].frontierSize; j++)
            G->inFrontier[G->TD[G->activeScatter[i]].frontier[j]] = false;
        G->TD[G->activeScatter[i]].frontierSize = 0;
    }
    G->frontierSize = 0;
}

template<class graph>
void loadFrontier(graph* G, intV* initFrontier, intV initFrontierSize)
{
    partitionData* allTD = G->TD;
    intV vertexId, pId;
    G->partListPtr = 0;
    for (intV i=0; i<initFrontierSize; i++)
    {
        vertexId = initFrontier[i];
        G->inFrontier[vertexId] = true;
        pId = (vertexId >> binOffsetBits);
        allTD[pId].frontier[allTD[pId].frontierSize++] = vertexId;
        allTD[pId].activeEdges += G->outDeg[vertexId];
        if (!G->flag[pId])
        {
            G->flag[pId] = true;
            G->activeScatter[G->partListPtr++] = pId;
        }
    }
    for (intV i=0; i<initFrontierSize; i++)
    {
        pId = (initFrontier[i] >> binOffsetBits);
        if (G->flag[pId])
            G->flag[pId] = false;
    }
    G->frontierSize = initFrontierSize;
}

/* 初始化图中的活跃节点信息
 * - 初始化全局活跃节点
 * - 初始化每个分块对应的活跃节点、活跃边数 
 * - 初始化分块的活跃情况和需要scatter的活跃块集合 */
template<class graph>
void loadFrontierPar(graph* G, intV* initFrontier, intV initFrontierSize)
{
    partitionData* allTD = G->TD;
    G->partListPtr = 0;
    #pragma omp parallel for
    for (intV i=0; i<initFrontierSize; i++)
    {
        intV vertexId = initFrontier[i];
        G->inFrontier[vertexId] = true;
        intV pId = (vertexId >> binOffsetBits);
        intV ptr = __sync_fetch_and_add(&allTD[pId].frontierSize, 1);
        allTD[pId].frontier[ptr] = vertexId;
        __sync_fetch_and_add(&allTD[pId].activeEdges, G->outDeg[vertexId]);
        if (__sync_bool_compare_and_swap(&G->flag[pId], false, true))
        {
            ptr = G->partListPtr.fetch_add(1);
            G->activeScatter[ptr] = pId;
        }
    }
    #pragma omp parallel for
    for (intV i=0; i<initFrontierSize; i++)
    {
        intV pId = (initFrontier[i] >> binOffsetBits);
        if (G->flag[pId])
            G->flag[pId] = false;
    }
    G->frontierSize = initFrontierSize;
}

/* 重新初始化sparse frontier数组, 对所有的活跃节点判断一轮 */
template <class graph, class userArg>
void reInitializeSparseFrontier(graph* G, partitionData* TD, userArg UA)
{
    intV trueSize = 0;
    for (intV i=0; i<TD->frontierSize; i++)
    {
        G->inFrontier[TD->frontier[i]] = UA.initFunc(TD->frontier[i]);
        if (G->inFrontier[TD->frontier[i]])
        {
            TD->frontier[trueSize++] = TD->frontier[i];
        }
    }
    if ((trueSize > 0) && 
        (__sync_bool_compare_and_swap(&G->flag[TD->tid], false, true)))
    {
       intV listPtr = G->partListPtr.fetch_add(1); 
       G->activeGather[listPtr] = TD->tid; 
    }
    TD->frontierSize = trueSize;
}

template <class graph, class userArg>
void reInitializeDenseFrontier(graph* G, partitionData* TD, userArg UA)
{
    for (intV i=TD->startVertex; i<TD->endVertex; i++)
    {
        UA.initFunc(i);
    }
}

/* gather之后重新构建全图的活跃节点信息和结构数组 */
template <class graph, class userArg>
void filterFrontier(graph* G, partitionData* TD, userArg UA)
{
    intV trueSize = 0;
    for (intV i=0; i<TD->frontierSize; i++)
    {
        G->inFrontier[TD->frontier[i]] = UA.filterFunc(TD->frontier[i]);
        if (G->inFrontier[TD->frontier[i]])
        {
            TD->activeEdges += G->outDeg[TD->frontier[i]];
            TD->frontier[trueSize++] = TD->frontier[i];
        }
    }
    TD->frontierSize = trueSize;
    /*
    if(G->numIter == 0){
        std::cout << "<<<<<<<<<<<<<<<<" << std::endl;
        std::cout << "frontier true size = " <<TD->frontierSize << std::endl;
        for(unsigned i = 0;i < TD->frontierSize;i++)
            std::cout << TD->frontier[i] << " ";
        std::cout << std::endl;
    }
    */
    if (TD->frontierSize > 0)
        G->frontierSize.fetch_add(TD->frontierSize);
}

/* 检查分块计算量，选择计算模式 */
void densityCheck(partitionData* TD)
{
    #ifdef SPARSE_MODE
        TD->isDense = false;
    #else
        TD->isDense = ((28.0 * (float)TD->activeEdges) > ((float)(TD->PNG->numEdges)*10.5 + 2.67*(float)TD->totalEdges + 4.0*(float)NUM_BINS)); 
    #endif
}

/* scatterVC函数:
 * - TD: 某一个part
 * - updateBins: 保存更新值的bin
 * - destIdBins: 保存目标节点id的bin
 * - updateBinPointers: 保存目前更新值的数量
 * - destIdBinPointers: 保存目标节点id值的数量
 * - UA: 用户定义的更新函数
 * - 这个函数在做正常的scatter操作, 这个过程中涉及bin的随机访问*/
template <class type,class graph, class userArg>
void scatterVC(graph* G, partitionData* TD, 
                type** updateBins, intV** destIdBins, 
                intE* updateBinPointers, intE* destIdBinPointers, userArg UA)
{
    intV destId = 0;
    intV destBin = 0;
    intV vertexId = 0;
    intV cond = 0;
    intV prevBin = 0;
    intV listPtr = 0;
    type userReturn;

#ifdef WEIGHTED
    type weightedVal;
#endif
    /* updateBin和destIdBin机制解读:
     * - updateBin中存放更新信息对应的更新值
     * - destBins中存放目标节点值, 目标节点为负的时候代表一个更新信息结束
     */
    for (intV i=0; i<TD->frontierSize; i++)
    {
        vertexId = TD->frontier[i]; //pop an active vertex
        prevBin = NUM_BINS;
        userReturn = UA.scatterFunc(vertexId); //invoke user def func. on the vertex 
        for (intE j=G->VI[vertexId]; j<G->VI[vertexId+1]; j++)
        {
            destId = G->EI[j]; 
            //if(vertexId == G->start)
            //    std::cout << destId << " " << std::endl;
            destBin = (destId >> binOffsetBits);

#ifdef WEIGHTED
            weightedVal = UA.applyWeight(userReturn, G->EW[j]); //apply weight to the value 
            updateBins[destBin][destIdBinPointers[destBin]] = weightedVal; //store the update in update bins
            destIdBins[destBin][destIdBinPointers[destBin]++] = destId; //store the dest ID in destination bins
#else
            ///////////////////////////////////
            ///// branch avoiding approach ////
            ///////////////////////////////////
            cond = (destBin != prevBin);
            updateBins[destBin][updateBinPointers[destBin]] = userReturn;
            updateBinPointers[destBin] += cond;
            destId |= (cond << MSB_ROT);
            destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
            prevBin = destBin;

            ///////////////////////////////////
            //////// branched approach ////////
            ///////////////////////////////////
            //if (destBin!=prevBin)
            //{
            //    updateBins[destBin][updateBinPointers[destBin]++] = userReturn;
            //    destId |= (1 << MSB_ROT);
            //    prevBin = destBin;
            //}
            //destIdBins[destBin][destIdBinPointers[destBin]++] = destId;

           
#endif
            if (!G->binFlag[TD->tid][destBin]) 
            {
                /* 如果对应的消息bin目前还没有被标记为活跃的话:
                 * - 标记该bin的状态为活跃
                 * - 为该目标part要gather的bin数量加一
                 * - 在二维的活跃bin列表中添加这个bin的ID
                 * - 如果对应的高层front还没有标记为活跃的话，则标记为活跃*/
                G->binFlag[TD->tid][destBin] = true; //mark the bin as active
                listPtr = G->TD[destBin].binListPtr.fetch_add(1); //atomically increase # active bins for destination partition
                G->activeBins[destBin][listPtr] = TD->tid; //convey the ID of bin to destination partition
                if (__sync_bool_compare_and_swap(&G->flag[destBin], false, true)) //if the destination partition isn't active
                {
                    listPtr = G->partListPtr.fetch_add(1); //
                    G->activeGather[listPtr] = destBin; 
                }
            }
        }
    }
}

template <class type,class graph, class userArg>
void scatterVCRW(graph* G, partitionData* TD, 
                type** updateBins, intV** destIdBins, 
                intE* updateBinPointers, intE* destIdBinPointers, userArg UA)
{
    intV destId = 0;
    intV destBin = 0;
    intV vertexId = 0;
    intV cond = 0;
    intV prevBin = 0;
    intV listPtr = 0;
    type userReturn;

#ifdef WEIGHTED
    type weightedVal;
#endif
    /* updateBin和destIdBin机制解读:
     * - updateBin中存放更新信息对应的更新值
     * - destBins中存放目标节点值, 目标节点为负的时候代表一个更新信息结束
     */
    for (intV i=0; i<TD->frontierSize; i++)
    {
        vertexId = TD->frontier[i]; //pop an active vertex
        prevBin = NUM_BINS;
        for(int num_msg = 0; num_msg < G->rw_prop[vertexId]; num_msg++){
            userReturn = UA.scatterFunc(vertexId); //invoke user def func. on the vertex 
            for (intE j=G->VI[vertexId]; j<G->VI[vertexId+1]; j++)
            {
                destId = G->EI[j]; 
                //if(vertexId == G->start)
                //    std::cout << destId << " " << std::endl;
                destBin = (destId >> binOffsetBits);

    #ifdef WEIGHTED
                weightedVal = UA.applyWeight(userReturn, G->EW[j]); //apply weight to the value 
                updateBins[destBin][destIdBinPointers[destBin]] = weightedVal; //store the update in update bins
                destIdBins[destBin][destIdBinPointers[destBin]++] = destId; //store the dest ID in destination bins
    #else
                ///////////////////////////////////
                ///// branch avoiding approach ////
                ///////////////////////////////////
                cond = (destBin != prevBin);
                updateBins[destBin][updateBinPointers[destBin]] = userReturn;
                updateBinPointers[destBin] += cond;
                destId |= (cond << MSB_ROT);
                destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
                prevBin = destBin;

                ///////////////////////////////////
                //////// branched approach ////////
                ///////////////////////////////////
                //if (destBin!=prevBin)
                //{
                //    updateBins[destBin][updateBinPointers[destBin]++] = userReturn;
                //    destId |= (1 << MSB_ROT);
                //    prevBin = destBin;
                //}
                //destIdBins[destBin][destIdBinPointers[destBin]++] = destId;

            
    #endif
                if (!G->binFlag[TD->tid][destBin]) 
                {
                    /* 如果对应的消息bin目前还没有被标记为活跃的话:
                    * - 标记该bin的状态为活跃
                    * - 为该目标part要gather的bin数量加一
                    * - 在二维的活跃bin列表中添加这个bin的ID
                    * - 如果对应的高层front还没有标记为活跃的话，则标记为活跃*/
                    G->binFlag[TD->tid][destBin] = true; //mark the bin as active
                    listPtr = G->TD[destBin].binListPtr.fetch_add(1); //atomically increase # active bins for destination partition
                    G->activeBins[destBin][listPtr] = TD->tid; //convey the ID of bin to destination partition
                    if (__sync_bool_compare_and_swap(&G->flag[destBin], false, true)) //if the destination partition isn't active
                    {
                        listPtr = G->partListPtr.fetch_add(1); //
                        G->activeGather[listPtr] = destBin; 
                    }
                }
            }
        }
    }
}


/* scatterPC函数:
 * - TD: 某一个part
 * - updateBins: 保存更新值的bin
 * - 这个函数的scatter操作是按照PNG存放的图格式写入更新*/
template <class type,class graph, class userArg>
void scatterPC(graph*G, partitionData* TD, type** updateBins, userArg UA)
{
    partitionGraph* PNG = TD->PNG;
    intE pointer;
    intV listPtr;
    type userReturn;
#ifndef DENSE
    for (intV i=0; i<NUM_BINS; i++)
    {
        G->binFlag[TD->tid][i] = true;
        listPtr = G->TD[i].binListPtr.fetch_add(1);
        G->activeBins[i][listPtr] = TD->tid;
        if (__sync_bool_compare_and_swap(&G->flag[i], false, true))
        {
            listPtr = G->partListPtr.fetch_add(1);
            G->activeGather[listPtr] = i; 
        }
    }
#endif
    for (intV i=0; i<PNG->numVertex; i++)
    {
        pointer = 0;
        for (intE j=PNG->VI[i]; j<PNG->VI[i+1]; j++){
            userReturn = UA.scatterFunc(PNG->EI[j]);
            updateBins[i][pointer++] = userReturn;
        }
    } 
}

/* scatter函数:
 * - TD: 某一个part
 * - updateBins: 保存更新的bin
 * - destIdBins: 保存目标节点id的bin
 * - updateBinPointers: 保存下一个bin地址
 * - destIdBinPointers: 保存下一个目标节点id地址
 * - UA: 用户定义的更新函数*/
template <class type, class graph, class userArg>
void scatter(graph*G, partitionData* TD, 
            type** updateBins, intV** destIdBins, 
            intE* updateBinPointers, intE* destIdBinPointers, userArg UA)
{
#ifndef DENSE
    if (TD->isDense){
        scatterPC<type>(G, TD, updateBins, UA);
    }
    else{
        scatterVC<type>(G, TD, updateBins, destIdBins, updateBinPointers, destIdBinPointers, UA);
        //scatterVCRW<type>(G, TD, updateBins, destIdBins, updateBinPointers, destIdBinPointers, UA);
    }
    reInitializeSparseFrontier(G, TD, UA);
#else
    scatterPC<type>(G, TD, updateBins, UA);
    reInitializeDenseFrontier(G,TD,UA);
#endif
}


//////////////////////////////////
/////////// PC GATHER ////////////
//////////////////////////////////
/* gatherPC参数:
 * - updateBin: 保存更新信息中更新的bin
 * - destIdBin: 保存更新信息中目标id的bin
 * - binSize: 更新条数
 * - 该函数的作用: 将更新写到目标节点上
 * */
template <class type, class graph, class userArg>
#ifdef WEIGHTED
void gatherPC(graph* G, partitionData* TD, type* updateBin, intV* destIdBin, unsigned int* weightBin, intE binSize, userArg UA)
#else
void gatherPC(graph* G, partitionData* TD, type* updateBin, 
                intV* destIdBin, intE binSize, userArg UA)
#endif
{
    intV destId = 0;
    intE updateBinPointer = MAX_UINT;
    type updateVal;
    bool cond;
    for (intE j=0; j<binSize; j++)
    {
        destId = destIdBin[j];
        // std::cout << std::endl;
        // 如果destId为负数的话, 指针向后移动1
        updateBinPointer += (destId >> MSB_ROT); 
        destId = destId & MAX_POS;
        
#ifdef WEIGHTED
        updateVal = UA.applyWeight(updateBin[updateBinPointer], weightBin[j]);
#else
        updateVal = updateBin[updateBinPointer];
#endif
        cond = UA.gatherFunc(updateVal, destId);
        // 更新frontier信息
        if (!G->inFrontier[destId] && cond)
        {
            TD->frontier[TD->frontierSize++] = destId;
            G->inFrontier[destId] = true;
        }
   } 
}

///////////////////////////////////
///////// DENSE GATHER ///////////
//////////////////////////////////
template <class type, class graph, class userArg>
#ifdef WEIGHTED
void gatherDense(graph* G, partitionData* TD, type* updateBin, intV* destIdBin, unsigned int* weightBin, intE binSize, userArg UA)
#else
void gatherDense(graph* G, partitionData* TD, type* updateBin, intV* destIdBin, intE binSize, userArg UA)
#endif
{
    intV destId = 0; 
    intE updateBinPointer = MAX_UINT;
    type updateVal;
    for (intE j=0; j<binSize; j++)
    {
        destId = destIdBin[j];
        updateBinPointer += (destId >> MSB_ROT);
        destId = destId & MAX_POS;
#ifdef WEIGHTED
        updateVal = UA.applyWeight(updateBin[updateBinPointer], weightBin[j]);
#else
        updateVal = updateBin[updateBinPointer];  
#endif
        UA.gatherFunc(updateVal, destId); 
    } 
}

#ifdef WEIGHTED
///////////////////////////////////
///////// VC GATHER ///////////
//////////////////////////////////
template <class type, class graph, class userArg>
void gatherVC(graph* G, partitionData* TD, type* updateBin, intV* destIdBin, intE binSize, userArg UA)
{
    
    intV destId = 0; 
    intE updateBinPointer = MAX_UINT;
    bool cond;
    type updateVal;
    for (intE j=0; j<binSize; j++)
    {
        destId = destIdBin[j];
        updateVal = updateBin[j];  
        cond = UA.gatherFunc(updateVal, destId);
        if (!G->inFrontier[destId] && cond)
        {
            TD->frontier[TD->frontierSize++] = destId;
            G->inFrontier[destId] = true;
        }
    } 
}
#endif


/* gather函数:
 * - updateBins: 更新信息的bin(对dense和sparse一致)
 * - denseDestIdBins: 目标id的bin (dense, fixed order)
 * - sparseDestIdBins: 稀疏的目标id的bin (empty, to be filled on scatter)
 * - allTD: 分块数据结构数组
 * - destIdBinAddrSize: bin内含有的目标节点数量的二维数组
 * - destIdBinPointers: bin内下一个目标节点写入地址的二维数组
 * - udpateBinPointers: bin内下一个更新写入地址的二维数组*/
template <class type, class graph, class userArg>
void gather(graph* G, partitionData* TD, 
            type*** updateBins, intV*** denseDestIdBins, 
            intV*** sparseDestIdBins, partitionData* allTD, 
            intE** destIdBinAddrSize, intE** destIdBinPointers, 
            intE** updateBinPointers, userArg UA)
{
    TD->activeEdges = 0;
#ifndef DENSE
    G->flag[TD->tid] = false;
#endif
#ifdef WEIGHTED
    unsigned int*** weightBin = G->indWeightBins;
#endif

    for (intV ptr=0; ptr<TD->binListPtr; ptr++)
    {
#ifndef DENSE
        intV i = G->activeBins[TD->tid][ptr];
        if (!G->binFlag[i][TD->tid]) continue; //already done during scatter-gather mix
        G->binFlag[i][TD->tid] = false;
#else
        intV i = ptr;
#endif
#ifdef WEIGHTED
#ifndef DENSE
        if (allTD[i].isDense)
            gatherPC<type>(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], weightBin[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
        else
            gatherVC<type>(G, TD, updateBins[i][TD->tid], sparseDestIdBins[i][TD->tid], destIdBinPointers[i][TD->tid], UA); 
#else
        gatherDense(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], weightBin[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
#endif
#else
#ifndef DENSE
        if (allTD[i].isDense)
            gatherPC<type>(G, TD, updateBins[i][TD->tid], 
                            denseDestIdBins[i][TD->tid], 
                            destIdBinAddrSize[i][TD->tid], 
                            UA);
        else
            gatherPC<type>(G, TD, updateBins[i][TD->tid], 
                            sparseDestIdBins[i][TD->tid], 
                            destIdBinPointers[i][TD->tid], 
                            UA); 
#else
        gatherDense(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
#endif
#endif
        destIdBinPointers[i][TD->tid] = 0;
        updateBinPointers[i][TD->tid] = 0;
    }
#ifndef DENSE
    TD->binListPtr = 0;
    filterFrontier(G, TD, UA);
#else
    for (intV i=TD->startVertex; i<TD->endVertex; i++)
        UA.filterFunc(i);
#endif
} 



template <class type, class graph, class userArg>
void gatherIL(graph* G, partitionData* TD, type*** updateBins, intV*** denseDestIdBins, intV*** sparseDestIdBins, partitionData* allTD, intE** destIdBinAddrSize, intE** destIdBinPointers, intE** updateBinPointers, bool* scatterDone, userArg UA)
{
    TD->activeEdges = 0;
#ifdef WEIGHTED
    unsigned int*** weightBin = G->indWeightBins;
#endif
    for (intV ptr=0; ptr<TD->binListPtr; ptr++)
    {
#ifndef DENSE
        intV i = G->activeBins[TD->tid][ptr];
#else
        intV i = ptr;
#endif
        if (!scatterDone[i]) continue;
        G->binFlag[i][TD->tid] = false;
#ifdef WEIGHTED
#ifndef DENSE
        if (allTD[i].isDense)
            gatherPC<type>(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], weightBin[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
        else
            gatherVC<type>(G, TD, updateBins[i][TD->tid], sparseDestIdBins[i][TD->tid], destIdBinPointers[i][TD->tid], UA); 
#else
        gatherDense(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], weightBin[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
#endif
#else
#ifndef DENSE
        if (allTD[i].isDense)
            gatherPC<type>(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
        else
            gatherPC<type>(G, TD, updateBins[i][TD->tid], sparseDestIdBins[i][TD->tid], destIdBinPointers[i][TD->tid], UA); 
#else
        gatherDense(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
#endif
#endif
        destIdBinPointers[i][TD->tid] = 0;
        updateBinPointers[i][TD->tid] = 0;
    }
#ifndef DENSE
    filterFrontier(G, TD, UA);
#else
    for (intV i=TD->startVertex; i<TD->endVertex; i++)
        UA.filterFunc(i);
#endif
} 

//for intra partition asynch processing
template <class type, class graph, class userArg>
void sgIntra(graph* G, partitionData* TD, userArg UA)
{
    partitionGraph* IGSort = TD->IPG;
    if (TD->isDense)
    {
        for (intV i=0; i<IGSort->numVertex; i++)
        {
            intV vertexId = TD->startVertex + i;
            type userReturn = UA.scatterFunc(vertexId);            
            for (intE j=IGSort->VI[i]; j<IGSort->VI[i+1]; j++)
            {
                intV destId = IGSort->EI[j];
#ifdef WEIGHTED
                type updateVal = UA.applyWeight(userReturn, IGSort->EW[j]);
#else
                type updateVal = userReturn;
#endif                 
                bool cond = UA.gatherFunc(updateVal, destId);
                if (!G->inFrontier[destId] && cond)
                {
                    TD->frontier[TD->frontierSize++] = destId;
                    G->inFrontier[destId] = true;
                }
            }
        }
    }
    else
    {
        for (intV i=0; i<TD->frontierSize; i++)
        {
            intV vertexId = TD->frontier[i];
            intV inPartId = vertexId - TD->startVertex;
            type userReturn = UA.scatterFunc(vertexId);
            for (intE j=IGSort->VI[inPartId]; j<IGSort->VI[inPartId+1]; j++)
            {
                intV destId = IGSort->EI[j];
#ifdef WEIGHTED
                type updateVal = UA.applyWeight(userReturn, IGSort->EW[j]);
#else
                type updateVal = userReturn;
#endif                 
                bool cond = UA.gatherFunc(updateVal, destId);
                if (!G->inFrontier[destId] && cond)
                {
                    TD->frontier[TD->frontierSize++] = destId;
                    G->inFrontier[destId] = true;
                }
            }
        }
    }
}


template <class type, class graph, class userArg>
void sgMix(graph* G, partitionData* TD, type*** updateBins, intV*** denseDestIdBins, intV*** sparseDestIdBins, partitionData* allTD, intE** destIdBinAddrSize, intE** destIdBinPointers, intE** updateBinPointers, bool* scatterDone, userArg UA)
{
    sgIntra<type>(G, TD, UA);
    gatherIL<type>(G, TD, updateBins, denseDestIdBins, sparseDestIdBins, allTD, destIdBinAddrSize, destIdBinPointers, updateBinPointers, scatterDone, UA); 
    densityCheck(TD);
    scatter<type>(G, TD, updateBins[TD->tid], sparseDestIdBins[TD->tid], updateBinPointers[TD->tid], destIdBinPointers[TD->tid], UA);
    scatterDone[TD->tid] = true;
}
