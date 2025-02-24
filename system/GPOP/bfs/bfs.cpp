/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * This code implements work efficient BFS
 * 
 */

unsigned int numIter = 0;

#include "../include/pcp.h"

struct BFS_F{
    intV* parent;
    bool* visited;
    BFS_F(intV* _parent, bool* _visited):parent(_parent), visited(_visited){}
    inline intV scatterFunc (intV node)
    {
        return (((!visited[node])<<MSB_ROT) | node);
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (intV updateVal, intV destId)
    {
        if (!visited[destId]) //if destination vertex is not yet visited
        {
            parent[destId] = updateVal; //set its parent
            visited[destId] = (!(updateVal>>MSB_ROT)); //new visited status depends on parent's status
            return visited[destId]; //active if it is now visited
        }
        return false;
    }  
    
    inline bool filterFunc(intV node)
    {
        return true;
    } 

};

struct BFS_F2{
    intV* parent;
    BFS_F2(intV* _parent):parent(_parent){}
    inline intV scatterFunc (intV node)
    {
        return ((parent[node] & MAX_NEG)| node);
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (intV updateVal, intV destId)
    { 
        //if parent of destID is negative and received update from a visited node
        if((parent[destId] & MAX_NEG) && (!(updateVal>>MSB_ROT)))
        {
            parent[destId] = updateVal; //update parent
            return true;
        }
        return false;
    }  
    
    inline bool filterFunc(intV node)
    {
        return true;
    } 
};


int main(int argc, char** argv)
{
    graph<intV> G;
    initialize(&G, argc, argv);
    initBin<intV>(&G);
    intV n = G.numVertex;
    intV* parent = new intV [n]();
    bool* visited = new bool [n]();
    struct timespec start, end;
    float time;

    int ctr = 0;
    float T = 0;
#ifdef SAMPLE_MODE
    printf("-----( Sampling Mode )-----\n");
    getRandStart(&G, G.rounds);
#endif

    while(ctr < G.rounds){
#ifdef SAMPLE_MODE
        std::string log_dir = G.dir + "/sample";
        // std::cout << "sampling dir = " << log_dir << std::endl;
        // std::cout << "data name = " << G.data_name << std::endl;
        if(!directoryExists(log_dir)){
            if(!createDirectory(log_dir)){
                std::cout << "Can not create directory " << log_dir << std::endl;
            }
        }
        std::string log_file = log_dir + "/" + G.data_name + ".bfs." + std::to_string(ctr);
        std::ifstream file(log_file.c_str());
        if(file.good()){
            std::remove(log_file.c_str());
        }
        file.close();
        std::ofstream ofs(log_file.c_str(), std::ios::app);
        ofs << G.numVertex << std::endl;
        G.start = G.starters[ctr];
        printf("[ round %d ]: starter = %d, sample file preserved at %s\n", ctr, G.start, log_file.c_str());
        ofs.close();
#endif
        intV initFrontierSize = 1;
        intV* initFrontier = new intV [initFrontierSize];
        for (intV i=0; i<initFrontierSize; i++)
            initFrontier[i] = G.start;
        for(int i=0;i<n;i++)
        {
            parent[i] = 1 << MSB_ROT;
            visited[i] = false;
        }
        visited[G.start] = true;
        parent[G.start] = G.start;
        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
    
        loadFrontierPar(&G, initFrontier, initFrontierSize);

        G.numIter = 0;
        G.scatterT = 0.0;
        G.gatherT = 0.0;
        unsigned int scatter_num = 0;

        while((G.frontierSize > 0))
        {
            scatter_num += G.partListPtr;
            #ifdef SAMPLE_MODE
                scatter_and_gather<intV>(&G, BFS_F(parent, visited), log_file);
            #else
                scatter_and_gather<intV>(&G, BFS_F(parent, visited), "");
            #endif

            G.numIter++;
        }
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("bfs, iter = %d, %lf\n", G.numIter, time);
        ctr++;
        T += time; 
    }
    printf("Avg time = %lf\n", T / G.rounds);
#ifdef SAMPLE_MODE
    printf("----------\nSample finished!\n");
#endif
    return 0;
}

    



