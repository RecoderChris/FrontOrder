#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <unordered_map>
#include "graph.h"

#define DEBUG 
#undef DEBUG

bool weighted = false;
bool createReverse = false;
bool undirected = false;

using namespace std;

int main(int argc, char** argv)
{

    unsigned int wtValue = 0;
    unsigned int reverseFileIndex = 0;

    if (argc < 3)
    {
        printf("Usage : %s <inputFile1> <outputFile> -w <weight_type> -r <graphTransposeOutputFilename>\n", argv[0]);
        exit(1);
    }
    // read csr
    graph G1;
    read_csr (argv[1], &G1);
    FILE* fp = fopen (argv[2], "w");
    for (intV src=0; src<G1.numVertex; src++){
        intE beg = G1.VI[src];
        intE end = G1.VI[src+1];
        for(intE idx=beg; idx < end;idx++){
            intV dst = G1.EI[idx];
            fprintf(fp, "%d %d\n", src, dst);
        }

    }
    fclose(fp);


    return 0;
}

