/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * Weakly connected components
 * 
 */

unsigned int numIter = 0;

//for asynchronous update propagation//
//converges faster//
//#define ASYNCH

#include "../include/pcp.h"

struct CC_F{
    intV* label;
    CC_F(intV* _label):label(_label){}

    inline intV scatterFunc (intV node)
    {
        return label[node];
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (intV updateVal, intV destId)
    {
        bool cond = (updateVal < label[destId]);
        if (cond)
            label[destId] = updateVal;
        return cond;
    }  
    
    inline bool filterFunc(intV node)
    {
        return true;
    } 

};


int main(int argc, char** argv)
{
    graph<intV> G;
    string data_name = initialize(&G, argc, argv);
    initBin<intV>(&G);
    intV n = G.numVertex;
    intV* label = new intV [n]();
    intV initFrontierSize = n;
    intV* initFrontier = new intV [initFrontierSize];
    for (intV i=0; i<initFrontierSize; i++)
        initFrontier[i] = i;

    struct timespec start, end;
    float time;
    
    /*
    while((G.frontierSize > 0))
    {
         #ifdef SAMPLE_MODE
            scatter_and_gather<intV>(&G, CC_F(label), log_file);
        #else
            scatter_and_gather<intV>(&G, CC_F(label), "");
        #endif
        numIter++;
    }
    */
    numIter = 0;
    int ctr = 0;
    float T = 0;
    while(ctr < G.rounds){
#ifdef SAMPLE_MODE
        std::string log_dir = "/home/zhangxm/project/data/" + data_name + "/sample/";
        std::cout << "sampling dir = " << log_dir << std::endl;
        if(!directoryExists(log_dir)){
            if(!createDirectory(log_dir)){
                std::cout << "Can not create directory " << log_dir << std::endl;
            }
        }
        std::string log_file = log_dir + data_name + ".cc."  + std::to_string(ctr);
        std::cout << "log file name = " << log_file << std::endl;
        std::ifstream file(log_file.c_str());
        if(file.good()){
            if(std::remove(log_file.c_str()) == 0)
                std::cout << "文件删除成功" << std::endl;
            else
                std::cout << "文件无法删除" << std::endl;
        }
        file.close();
        std::ofstream ofs(log_file.c_str(), std::ios::app);
        ofs << G.numVertex << std::endl;
        std::cout << "========== [round " << ctr << " ] ==========" << std::endl;
#endif
        for(int i=0;i<n;i++){
            label[G.new_id[i]] = i;
        }
        loadFrontier(&G, initFrontier, initFrontierSize);

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}


        while((G.frontierSize > 0))
        {
            // printf("-----[ Iter = %u ]-----\n", numIter);
            // printf("Frontier size = %u\n", G.frontierSize.load());
            #ifdef SAMPLE_MODE
                scatter_and_gather<intV>(&G, CC_F(label), log_file);
            #else
                scatter_and_gather<intV>(&G, CC_F(label), "");
            #endif
            numIter++;
        }
#ifdef SAMPLE_MODE
        ofs.close();
#endif
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        T += time;
        printf("cc, %d, %s, %lf\n", NUM_THREADS, argv[1], time);
        ctr++;
    }
    printf("Avg time = %lf\n", T / G.rounds);


    printf("\n");


    return 0;
}



