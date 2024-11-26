// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "ligra.h"
#include <numeric>
#define WLB
//#undef WLB
#define THREAD_NUM 96

std::vector<unsigned> workload(THREAD_NUM);


struct BFS_F {
  uintE* Parents;
  BFS_F(uintE* _Parents) : Parents(_Parents) {}
  inline bool update (uintE s, uintE d) { //Update
    
    if(Parents[d] == UINT_E_MAX) {  
      #ifdef OPENMP
        workload[omp_get_thread_num()]++;
      #endif
      Parents[d] = s; return 1; 
    }
    else return 0;
  }
  inline bool updateAtomic (uintE s, uintE d){ //atomic version of Update
    #ifdef OPENMP
      workload[omp_get_thread_num()]++;
    #endif
    return (CAS(&Parents[d],UINT_E_MAX,s));
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return (Parents[d] == UINT_E_MAX); } 
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long start = P.getOptionLongValue("-r",0);
  // std::cout << "orignal id = " << start << std::endl;
  start = map[start];
  // std::cout << "new id = " << map[10] << std::endl;
  // std::cout << "new id = " << start << std::endl;
  // assert(0);
  
  long n = GA.n;
  //creates Parents array, initialized to all -1, except for start
  uintE* Parents = newA(uintE,n);
  parallel_for(long i=0;i<n;i++) Parents[i] = UINT_E_MAX;
  Parents[start] = start;
  vertexSubset Frontier(n,start); //creates initial frontier
  // std::vector<double> wlb_factor = {0.0};
  // std::vector<unsigned> max_work;
  // std::vector<unsigned> avg_work;
  while(!Frontier.isEmpty()){ //loop until frontier is empty
    // std::fill(workload.begin(), workload.end(), 0);
    vertexSubset output = edgeMap(GA, Frontier, BFS_F(Parents));    
    Frontier.del();
    Frontier = output; //set new frontier
    // std::sort(workload.begin(), workload.end());
    // unsigned max_val = workload[workload.size() - 1];
    // max_work.push_back(max_val);
    // unsigned avg_val = std::ceil((double)std::accumulate(workload.begin(), workload.end(), 0) / workload.size());
    // avg_work.push_back(avg_val);
    // if(avg_val > 0)
        // wlb_factor.push_back((double)max_val / avg_val);
  }
  // double avg_wlb = (double)std::accumulate(wlb_factor.begin(), wlb_factor.end(), 0) / wlb_factor.size();
  // printf("[avg wlb] = %.4lf\n", avg_wlb);
  // for(int i = 0;i < max_work.size();i++){
  //   printf("\t[%d]: max = %d, avg = %d\n", i, max_work[i], avg_work[i]);
  // }
  // int sum_max = std::accumulate(max_work.begin(), max_work.end(), 0);
  // printf("[sum max] = %d\n", sum_max);

  Frontier.del();
  free(Parents); 
}
