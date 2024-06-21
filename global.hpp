#pragma once

#include <cmath>

//#define WEIGHTED

namespace params {
    float damping = 0.15;
    unsigned int partition_size = (1024 * 1024)/sizeof(float);//(256*1024)/sizeof(float); //512kB cluster size is for cluster constructing
    unsigned int num_partitions = 0;
    unsigned int partition_offset = (unsigned)log2((float)partition_size); 
    unsigned int num_threads = 10;
    unsigned int overflow_ceil = 0;
};


enum Algo {
    original = 0,
    hisorder_wo_blc = 1, 
    hisorder = 2, 
    Hisorder_cc_sc = 3, 
    HisOrder_cc = 4, 
    HisOrder_cc_noblc = 5, 
    HON_bfs = 6, 
    HO_bfs = 7, 
    bfs_blc = 8, 
    HO_SC = 12, 
    /*
    hisorder_it = 14,
    HO_bfs_thread = 7, 
    HO_bfs_mode = 8, 
    HON_cold = 9, 
    HO_cold = 10, 
    HO_mode = 11, 
    */
};