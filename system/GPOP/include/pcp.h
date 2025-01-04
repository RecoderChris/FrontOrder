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

#define CACHE_SIZE_KB  512

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <immintrin.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../include/gas.h"

#define CREATE_DIR(path) mkdir(path, 0777)
// #define SAMPLE_MODE
// #undef SAMPLE_MODE


#define PRINT_DENSE
#undef PRINT_DENSE

#define PRINT_WLC
#undef PRINT_WLC

intV binWidth = (CACHE_SIZE_KB * 1024)/sizeof(float); //512kB
unsigned int binOffsetBits = (unsigned int)std::log2((double)binWidth); 
intV NUM_BINS = 10000000/binWidth;


#define DEBUG
#undef DEBUG

//////////////////////////////////////////
// performance monitoring via PCM
//////////////////////////////////////////
#define PERF_MON
#undef PERF_MON


//////////////////////////////////////////
// level 2 debugging - asserts enabled
//////////////////////////////////////////
#define DEBUGL2
#undef DEBUGL2


#define ITERTIME
#undef ITERTIME

int NUM_THREADS = std::max(omp_get_max_threads(), 1);
unsigned int MAX_ITER = 100000;

bool directoryExists(const std::string& path) {
    struct stat info;
    return (stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode));
}

bool createDirectory(const std::string& path) {
    int result = CREATE_DIR(path.c_str());
    return (result == 0);
}

template<class graph>
void getRandStart(graph* G, unsigned r){
    unsigned v;
    srand((unsigned)time(NULL));
    for(int i = 0;i < r;i++){
        v = rand() % G->numVertex;
        while(G->outDeg[G->new_id[v]] < 1)
            v = rand() % G->numVertex;
        G->starters.push_back(v);
    }
}

template<class graph>
void initOrigOrder(graph* G) {
    G->new_id = std::vector<intV>(G->numVertex, 0);
    for (unsigned int i = 0 ; i < G->numVertex ; i++ )
        G->new_id[i] = i;
}

template<class graph>
void mapReorder(graph* G, std::string mapping_file) {
    std::ifstream ifs(mapping_file.c_str(), std::ifstream::in);
    if (!ifs.good()) {
        std::cout << "File " << mapping_file << " does not exist!" << std::endl;
        exit(-1);
    }
    unsigned long int num_vertex, num_edges;
    ifs >> num_vertex;
    ifs >> num_edges;
    G->new_id = std::vector<intV>(num_vertex, 0);
    // std::cout << " num_vertex: " << num_vertex << " num_edges: "  << num_edges << std::endl;
    char c;
    unsigned long int st, v, d2s;
    for ( unsigned int i = 0 ; i < num_vertex ; i++ ) {
        ifs >> st >> v;
        G->new_id[st] = v;
    }
    ifs.close();
}

template<class graph>
void initialize(graph* G, int argc, char** argv)
{
    string map_file = "";
    G->start = 1; 
    G->rounds= 5; 
    // std::cout << "=====[ Check Configs ]=====" << std::endl;
    for (int i = 1; i < argc; i++)
    {  
        if (i + 1 != argc)
        {
            if (strcmp(argv[i], "-s") == 0) // start node 
            {                 
                G->start = (intV)atoi(argv[i + 1]);    
                i++;
            }
            if (strcmp(argv[i], "-iter") == 0)  // num iterations
            {                 
                MAX_ITER = (unsigned int)atoi(argv[i + 1]);    
                i++;    
            }
            if (strcmp(argv[i], "-rounds") == 0) 
            {                 
                G->rounds = (unsigned int)atoi(argv[i + 1]);
                i++;   
            }
            if (strcmp(argv[i], "-map") == 0)
            {
                map_file = argv[i + 1];
                i++;
                // std::cout << "map file = " << map_file << std::endl;
                mapReorder(G, map_file);
            }
        }
    }
    
    if (argc < 2)
    {
        printf("Usage : %s <filename> -s <start node(not needed for pr)> -t <numThreads(optional)> -iter <#iterations(optional) -rounds <#rounds(default 3)> -map <map_file>\n", argv[0]);
        exit(1);
    }

    NUM_THREADS = omp_get_max_threads() / 2;
    printf("num threads = %d\n", NUM_THREADS);
    omp_set_num_threads(NUM_THREADS);

    //////////////////////////////////////////
    // read csr file
    //////////////////////////////////////////
    if (read_csr(argv[1], G)==-1)
    {
        printf("couldn't read %s\n", argv[1]);
        exit(1);
    }

    // Original graph order
    if(map_file == ""){
        initOrigOrder(G);
    }
#ifndef SAMPLE_MODE
    G->start = G->new_id[G->start];
    std::cout << "New starter = " << G->start << std::endl;
#endif // !
    // printf("L2 Cache(KB) = %d\n", CACHE_SIZE_KB);

    // intV numVerticesPerBin= binWidth;
    intV numVerticesPerBin= (G->numVertex/(NUM_THREADS*4));
    numVerticesPerBin = (numVerticesPerBin < binWidth) ? numVerticesPerBin : binWidth;
    // printf("CPU number = %d\n", NUM_THREADS);
    intV pow2=1;
    while(pow2<=numVerticesPerBin)
        pow2*=2;
    pow2/=2;
    if(pow2==0) binWidth=4;
    else binWidth = pow2;
    NUM_BINS = (G->numVertex-1)/binWidth + 1;
    G->numBins = NUM_BINS;
    binOffsetBits = (unsigned int)std::log2((double)binWidth);

    //////////////////////////////////////////
    //initialize graph frontier, degree etc.//
    //////////////////////////////////////////
    initGraph (G);
    // printf("Number of rounds (samplings) = %u\n", G->rounds);
}

template<class type, class graph>
void initBin(graph* G)
{
    //////////////////////////////////////////
    //static work allocation to threads
    //equal no. of edges to all bins
    //////////////////////////////////////////
    G->TD = (partitionData*) malloc (sizeof(partitionData)*NUM_BINS);
    partition(G->TD, G);

    //////////////////////////////////////////////////
    //compute storage space required for each bin and
    //offsets for storage in bins for a partition
    //1 column -> 1 gather bin; 1 row -> 1 scatter bin
    //bin[i][j] -> stores what i sends to j
    //////////////////////////////////////////////////
    G->updateBinAddrSize = allocateBinMat<intE>(NUM_BINS, NUM_BINS);
    G->destIdBinAddrSize = allocateBinMat<intE>(NUM_BINS, NUM_BINS);
    G->binFlag = allocateBinMat<bool>(NUM_BINS, NUM_BINS);
    G->activeBins = allocateBinMat<intV>(NUM_BINS, NUM_BINS);

//    struct timespec preStart, preEnd; 
//    float preTime;
//    if( clock_gettime(CLOCK_REALTIME, &preStart) == -1) { perror("clock gettime");}

    //////////////////////////////////////////
    //// transpose and compute offsets ///////
    //////////////////////////////////////////
    #pragma omp parallel for schedule (dynamic, 1)
    for (intV i=0; i<NUM_BINS; i++)
        transposePartition(G, &(G->TD[i]), G->updateBinAddrSize[i], G->destIdBinAddrSize[i]);




//    if( clock_gettime( CLOCK_REALTIME, &preEnd) == -1 ) { perror("clock gettime");}      
//    preTime = (preEnd.tv_sec - preStart.tv_sec)+ (int)(preEnd.tv_nsec - preStart.tv_nsec)/1e9;
//    printf("%s, preprocessing time - %lf\n", argv[1], preTime);

//////////////////////////////////////////
//////////////// BINNING ////////////////
//////////////////////////////////////////

    //////////////////////////////////////////
    ////individual bins to->fro each partition //////
    //////////////////////////////////////////
    G->indUpdateBins = allocateBinMatPtr<type>(NUM_BINS, NUM_BINS);
    G->indDestIdBins = allocateBinMatPtr<intV>(NUM_BINS, NUM_BINS);
    G->sparseDestIdBins = allocateBinMatPtr<intV>(NUM_BINS, NUM_BINS);
#ifdef WEIGHTED
    G->indWeightBins = allocateBinMatPtr<unsigned int>(NUM_BINS, NUM_BINS);
#endif
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, 1)
    for (intV i=0; i<NUM_BINS; i++)
    {
        for (intV j=0; j<NUM_BINS; j++)
        {
            G->indUpdateBins[i][j] = new type [G->destIdBinAddrSize[i][j]];
            G->indDestIdBins[i][j] = new intV [G->destIdBinAddrSize[i][j]];
            G->sparseDestIdBins[i][j] = new intV [G->destIdBinAddrSize[i][j]];
#ifdef WEIGHTED
            G->indWeightBins[i][j] = new unsigned int [G->destIdBinAddrSize[i][j]];
#endif
        }
    }

    //pointers for each (i,j) bin for later use //
    G->updateBinPointers = allocateBinMat<intE>(NUM_BINS, NUM_BINS);
    G->destIdBinPointers = allocateBinMat<intE>(NUM_BINS, NUM_BINS);


    #pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic, 1)
    for (intV i=0; i<NUM_BINS; i++)
    {
#ifdef WEIGHTED
        writeDestIds(G, &G->TD[i], G->indDestIdBins[i], G->indWeightBins[i], G->destIdBinPointers[i]);
#else
        writeDestIds(G, &G->TD[i], G->indDestIdBins[i], G->destIdBinPointers[i]);
#endif
    }

#ifdef DEBUG
    printf("binning complete\n");
#endif 

//////////////////////////////////////////
//////////// BINNING COMPLETE ////////////
//////////////////////////////////////////

}


template<class type, class graph, class userArg>
void scatter_and_gather(graph* G, userArg UA, std::string log_file)
{
    intV numActiveBins;
///////////////////////////////////////
////Set FLAG For Scatter and Gather////
///////////////////////////////////////


#ifndef DENSE
        G->frontierSize = 0;
#endif 

        numActiveBins = G->partListPtr;
#ifdef SAMPLE_MODE 
        std::ofstream ofs(log_file.c_str(), std::ios::app);
        if (!ofs) {
            std::cout << "Can not open the feature file" << std::endl;
            exit(-1);
        }
        ofs << "----------" << std::endl;
        for (intV i = 0; i < numActiveBins; i++) {
            partitionData* TD = &G->TD[G->activeScatter[i]];
            for(intV j = 0;j < TD->frontierSize;j++)
                ofs << TD->frontier[j] << " ";
        }
        ofs << std::endl;
        ofs.close();
#endif

#ifndef DENSE
        G->partListPtr = 0;
#endif
#ifndef ASYNCH
        #pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic, 1)
        for (intV i=0; i < numActiveBins; i++) {
            densityCheck(&G->TD[G->activeScatter[i]]);
        }
#endif
        #pragma omp parallel for schedule(dynamic,1) num_threads(NUM_THREADS)
        for (intV ptr=0; ptr < numActiveBins; ptr++)
        {
            intV i = G->activeScatter[ptr];
#ifdef ASYNCH
            sgMix(G, &G->TD[i], G->indUpdateBins, G->indDestIdBins, G->sparseDestIdBins, G->TD, G->destIdBinAddrSize, G->destIdBinPointers, G->updateBinPointers, G->scatterDone, UA);
#else
            scatter<type>(G, &G->TD[i], G->indUpdateBins[i], 
                            G->sparseDestIdBins[i], G->updateBinPointers[i], 
                            G->destIdBinPointers[i], UA);
#endif
        } 
        #pragma omp parallel for
        for (intV ptr=0; ptr<numActiveBins; ptr++)
            G->scatterDone[G->activeScatter[ptr]] = false; //reset scatter done status of partitions
        #pragma omp parallel for schedule(dynamic,1) num_threads(NUM_THREADS)
        for (intV ptr=0; ptr<G->partListPtr; ptr++)
        {
            intV i=G->activeGather[ptr];
            gather<type>(G, &G->TD[i], G->indUpdateBins, G->indDestIdBins, 
                         G->sparseDestIdBins, G->TD, G->destIdBinAddrSize, 
                         G->destIdBinPointers, G->updateBinPointers, UA);
            G->activeScatter[ptr] = i;
        }    
}



