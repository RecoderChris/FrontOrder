#pragma once

#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include <vector>
#include <cmath>
#ifndef SORT_HEADER
#include "sort.hpp"
#include "global.hpp"
#endif

#define DEBUG_ATTR 
//#undef  DEBUG_ATTR

#define DEBUG_DETAIL
#undef DEBUG_DETAIL


class Graph{
public:
    unsigned num_vertex; // 图节点数
    unsigned num_edges; // 图边数
    unsigned cluster_num; // 簇的数量
    unsigned num_vertex_deg_0;
    unsigned num_vertex_deg_1;
    std::vector<unsigned> row_index; // csr中row数组
    std::vector<unsigned> col_index; // csr中col数组
    std::vector<unsigned> out_degree; // 图的出度数组
    std::vector<int> in_degree; // 图的入度数组
    unsigned feat; // 图的特征维度
    std::vector<std::vector<unsigned>> attr; // 节点property数组
    std::string in_feat;
#ifdef WEIGHTED
    std::vector<unsigned> edge_weight; // 边权重数组
#endif
    Graph(unsigned num_vertex = 0, unsigned num_edges = 0):
            num_vertex(num_vertex), num_edges(num_edges){
                row_index.reserve(num_vertex + 1);
                col_index.reserve(num_edges);
            }
    ~Graph(){}
    
    /* 通过csr中row数组计算图中所有节点的出度 */
    void computeOutDegree() {
        out_degree = std::vector<unsigned>(num_vertex, 0);
        in_degree = std::vector<int>(num_vertex, 0);
        uint64_t total = 0;
        float avr = num_edges / num_vertex;
        
        #pragma omp parallel for schedule(static, 1024 * 256) num_threads(params::num_threads)
        for(unsigned i = 0; i < num_vertex; i++) {
            out_degree[i] = row_index[i + 1] - row_index[i];
        }
        num_vertex_deg_0 = num_vertex;
        num_vertex_deg_1 = 0;
        for(int i = 0;i < num_vertex;i++){
            for(unsigned j = row_index[i];j < row_index[i + 1];j++){
                unsigned v = col_index[j];
                in_degree[v]++;
                /*
                if(in_degree[v] == 0){
                    num_vertex_deg_0--;
                    num_vertex_deg_1++;
                    in_degree[v] = i + 1; // 记录deg=1节点的in-degree
                    continue;
                }
                else if(in_degree[v] > 0){
                    num_vertex_deg_1--;
                    in_degree[v] = -1;
                }
                */
            }
        }
        //std::cout << "[Value Check]: " << std::endl;
        //std::cout << "\t#v(deg = 0) = " << num_vertex_deg_0 << std::endl;
        //std::cout << "\t#v(deg = 1) = " << num_vertex_deg_1 << std::endl;
    }
   
    /* 将全图以row/col/weight的形式打印出来 */
    void printGraph(bool all = false) {
        std::cout << "num vertex: " << num_vertex << std::endl;
        std::cout << "num edges: " << num_edges <<std::endl;
        #ifdef WEIGHTED
        std::cout << "weighted graph" <<std::endl;
        #endif
        if(all) {
            for(auto it = row_index.begin(); it != row_index.end(); ++it) 
                std::cout << *it << " ";
            std::cout << std::endl;
            for(auto it = col_index.begin(); it != col_index.end(); ++it) 
                std::cout << *it << " ";
            std::cout << std::endl;
        }
    }

    /* 初始化图节点特征数组 */
    void initAttribute(unsigned int feat)
    {
        this->feat = feat;
        // 初始化attr特征数组, 尺寸为 vertex_num x feat大小
        attr = std::vector<std::vector<unsigned>>(num_vertex);
        for(unsigned int i = 0;i < num_vertex;i++)
            attr[i] = std::vector<unsigned>(feat, std::numeric_limits<unsigned int>::max());
        
        srand((unsigned)time(NULL));
        for(unsigned int i = 0; i < feat;i++){
            unsigned start = std::rand() % num_vertex; // 随机产生度数不为0的节点
            while(row_index[start + 1] - row_index[start] < 1)
                start = std::rand() % num_vertex; // 随机产生一个节点
            printf("rand start vertex(%d) = %d\n", i, start);
            std::vector<unsigned> front = {start};
            unsigned int level = 0;
            attr[start][i] = level;
            while(front.size() > 0){
                // std::cout << "level = " << level << ", front size = " << front.size() << "" << "\n";
                std::vector<unsigned> next_front = {};
                level++;
                for(unsigned int v: front){
                    for(unsigned int col_idx = row_index[v]; col_idx < row_index[v + 1];col_idx++){
                        unsigned int dst = col_index[col_idx];
                        if(attr[dst][i] > level){
                            attr[dst][i] = level;
                            next_front.push_back(dst);
                        }
                    }
                }
                front = next_front;
            }
        }
        return;
    }

    /* 初始化图节点特征数组 */
    void initAttributeFile(std::string feat_file, int feat_size)
    {
        std::cout << "feat file = " << feat_file << std::endl;
        std::ifstream file(feat_file); // 更改为实际的文件名
        if (!file.is_open()) {
            std::cout << "无法打开文件,程序退出" << std::endl;
            exit(1);
        }
        // this->cluster_num = 0;
        std::string line;
        while (std::getline(file, line)) {
            std::vector<unsigned> row; // 用于存储每行的数据
            std::stringstream ss(line);
            std::string token;
            int len = 0;
            while (std::getline(ss, token, ',') && len < feat_size) {
                unsigned value = (unsigned)(std::stoi(token));
                // if(value > this->cluster_num)
                    // this->cluster_num = value;
                row.push_back(value);
                len++;
            }
            attr.push_back(row);
            if(attr.size() == 1){
                std::cout << "feat size = " << row.size() << std::endl;
                this->feat = row.size();
            }
        }

        file.close();

        this->feat = feat;
        return;
    }

    /* bfs */
    void bfs(unsigned starter)
    {
        std::vector<unsigned> prop(num_vertex, std::numeric_limits<unsigned>::max());
        std::vector<unsigned> front = {starter};
        unsigned int level = 0;
        unsigned int block_id = 0;
        unsigned int block_start, block_end;
        unsigned int num_total = 0;
        while(front.size() > 0){
            block_id = 0;
            block_start = front[0];
            //printf("----------\n");
            std::sort(front.begin(), front.end());
            std::vector<unsigned> next_front = {};
            level++;
            for(unsigned int i = 0;i < front.size();i++){
                unsigned int v = front[i];
                if(i > 0 && (v != front[i - 1] + 1)){
                    block_end = front[i - 1];
                #ifdef DEBUG_DETAIL
                    printf("(%d):  %u - %u\n", block_id, block_start, block_end);
                #endif
                    block_id++;
                    block_start = v;
                }
                if(i == front.size() - 1){
                    block_id++;
                #ifdef DEBUG_DETAIL
                    printf("(%d):  %u - %u\n", block_id, block_start, v);
                #endif
                }
                
                for(unsigned int col_idx = row_index[v]; col_idx < row_index[v + 1];col_idx++){
                    unsigned int dst = col_index[col_idx];
                    if(prop[dst] > level){
                        prop[dst] = level;
                        next_front.push_back(dst);
                    }
                }
            }
            //printf("本轮迭代活跃节点数量 = %d\n", front.size());
            //printf("本轮迭代连续活跃块数量 = %d\n", block_id + 1);
            num_total += (block_id + 1);

            front = next_front;
        }
        printf("==========\n总的活跃块数量 = %d\n", num_total);
        return;
    }
};

