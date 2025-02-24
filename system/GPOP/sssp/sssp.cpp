/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * Single source shortest path using Bellman-Ford 
 *
 */

unsigned int numIter = 0;

#define WEIGHTED

//for asynchronous update propagation//
//converges faster//
#define ASYNCH
#undef ASYNCH

#include "../include/pcp.h"


struct SSSP_F{
    unsigned int* distance;
    SSSP_F(unsigned int* _distance):distance(_distance){}

    inline unsigned int scatterFunc (intV node)
    {
        return distance[node];
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (unsigned int updateVal, intV destId)
    {
        if(updateVal < distance[destId])
        {
            distance[destId] = updateVal;
            return true;
        }
        else
            return false;
    }  
    
    inline bool filterFunc(intV node)
    {
        return true;
    } 

    inline unsigned int applyWeight (unsigned int updateVal, unsigned int weight)
    {
        return updateVal + weight;
    }

};




int main(int argc, char** argv)
{
    graph<unsigned int> G;
    initialize(&G, argc, argv);
    initBin<unsigned int>(&G);    
    intV n = G.numVertex;
    unsigned int* distance = new unsigned int [n]();
    
    struct timespec start, end;
    float time;

    int ctr = 0;
    float T = 0;
    // std::cout << "Old start id = " << G.start << std::endl;
    // G.start = G.new_id[G.start];
    // std::cout << "Old start id = " << G.start << std::endl;

    // for(unsigned st = 0;st < NUM_STARTERS;st++){
    //     std::cout << "-----------------" << std::endl;
    //     std::cout << "Start = " << st << std::endl;
    //     G.start = G.new_id[G.starters[st]];
    //     std::cout << "New start id = " << G.start << std::endl;
    //     ctr = 0;
    //     T = 0;
    while(ctr < G.rounds){
        intV initFrontierSize = 1;
        intV* initFrontier = new intV [initFrontierSize];
        for (intV i=0; i<initFrontierSize; i++)
            initFrontier[i] = G.start;

        loadFrontier(&G, initFrontier, initFrontierSize);
        
        numIter=0;
        for(int i=0;i<n;i++)
            distance[i] = 1<<31;
        distance[G.start] = 0;

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while((G.frontierSize > 0) && (numIter < G.numVertex))
        {
            scatter_and_gather<unsigned int>(&G, SSSP_F(distance), "");
            numIter++;
        }

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("sssp, %d, %s, %lf\n", NUM_THREADS, argv[1], time);
        ctr++;
        T += time;
    }
    printf("Avg time = %lf\n", T / G.rounds);
    printf("\n");
    // }
    
    /*
    FILE* fp = fopen("dump.bin", "wb");
    fwrite(distance, sizeof(unsigned int), n, fp);
    fclose(fp);
    */
    return 0;
}



