#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <limits.h>
#include <parallel/algorithm>
#include <parallel/numeric>
#include <algorithm>
#include <random>
#include <cmath>

#include "graph.hpp"
#include "vec2d.hpp"
#include "global.hpp"
#include <boost/timer/timer.hpp>
#include <vector>
#define DEBUG
#undef DEBUG
using namespace boost::timer;

#define KMEANS_ITER  100

#define DEBUG_HORDER
#undef DEBUG_HORDER

// 节点类型
class Node {
public:
    int id;
    int cluster_id;
    unsigned src_id;
    std::vector<unsigned> feat;

    // 计算两个节点之间的距离
    unsigned calculate_diff(Node node) {
        unsigned diff = 0;
        for (size_t i = 0; i < feat.size(); i++){
            int m = (int)feat[i] - (int)node.feat[i];
            diff += std::abs(m);
            // diff += (feat[i] != node.feat[i]);
        }
        return diff;
    }

    unsigned calculate_diff_kmp(Node node){
        unsigned diff = 0;
        for (size_t i = 0; i < feat.size(); i++){
            long int m = (long int)feat[i] - (long int)node.feat[i];
            diff += std::abs(m);
            //diff += (feat[i] != node.feat[i]);
        }
        return diff;
    }

    unsigned calculate_diff_onehot(Node node) {
        unsigned diff = 0;
        for (size_t i = 0; i < feat.size(); i++){
            if(feat[i] == 0 || node.feat[i] == 0)
                diff += 1;
        }
        return diff;
    }

    // 计算节点和任意一个向量的距离
    unsigned calculate_diff_vec(std::vector<unsigned> feat_vec){
        unsigned diff = 0;
        for(size_t i = 0; i < feat.size();i++){
            int m = (int)feat[i] - (int)feat_vec[i];
            diff += std::abs(m);
    
        }
            // diff += (feat[i] != feat_vec[i]);
        return diff;
    }

};

// 定义比较函数，用于指定排序依据
bool compareByValue(const Node& a, const Node& b) {
    return a.cluster_id < b.cluster_id;
}

// 重排器数据结构
class Orderer {
   // typedef std::vector<std::vector<unsigned>> Vector2D;
    typedef std::pair<unsigned, unsigned> degree_id_pair; // 度数-节点id对
    unsigned num_vertex; // 节点数量
    unsigned num_edges; // 边数量
    unsigned num_partitions; // 分块数量
    unsigned num_levels; // 层次分块中层次的数量
    unsigned average_degree; // 平均度数
    std::vector<unsigned> levels; // 平均度数的等比数列, 用于某些算法进行判定

    Graph* graph;
public:
    std::vector<unsigned> new_id;
    std::vector<unsigned> new_dst_id;

    /* 初始化重排结构 */
    Orderer(Graph* g) {
        this->graph = g;
        num_vertex = graph->num_vertex;
        num_edges = graph->num_edges;
        num_partitions = params::num_partitions;
        new_id = std::vector<unsigned>(num_vertex, 0);
        new_dst_id = std::vector<unsigned>(num_vertex, 0);
        average_degree = num_edges/num_vertex;

        // level: aver_degree/2, aver_degree, aver_degree*2
        num_levels = (unsigned) log2((float) num_partitions) + 2;

        for(int i = 0; i < num_levels; i++)
            levels.push_back(average_degree * pow(2, i - 1));

        levels.back() = UINT_MAX;
    }

    // 初始版本的horder
    void HisOrder_wo_blc() {
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }
        std::vector<Node> centroids; // 质心数组
        unsigned num_clusters = this->graph->cluster_num;
        printf("分类数量 = %d\n", num_clusters);
        std::vector<std::vector<Node>> clusters(num_clusters); // 聚类数组

        /* KMeans++初始化 */
        int first_center_id = rand() % nodes.size();
        centroids.push_back(nodes[first_center_id]);
        for(unsigned i = 1; i < num_clusters;i++){
            if(i % 10  == 0)
                printf("Centroid %d init..\n", i);
            unsigned total_distance_sq = 0;
            std::vector<unsigned> distances(num_vertex, std::numeric_limits<unsigned>::max());
            
            #pragma omp parallel for
            for(unsigned j = 0; j < nodes.size();j++){
                for(Node& centroid : centroids){
                    unsigned distance = centroid.calculate_diff_kmp(nodes[j]);
                    distances[j] = std::min(distances[j], distance);
                }
                unsigned sq = distances[j] * distances[j];
                #pragma omp atomic
                    total_distance_sq += sq;
            }
            double randValue = std::rand() / (RAND_MAX + 1.0) * total_distance_sq;
            for(unsigned j = 0; j < num_vertex;j++){
                randValue -= distances[j] * distances[j];
                if(randValue <= 0){
                    centroids.push_back(nodes[j]);
                    break;
                }
            }
        }
        printf("---- [ 排序之前 ] ----:\n");
        for(unsigned i = 0;i < centroids.size();i++){
            printf("Center(%d): ", i);
            for(unsigned j = 0;j < centroids[i].feat.size();j++){
                printf("%u ", centroids[i].feat[j]);
            }
            printf("\n");
        }
        for(unsigned i = 0;i < centroids.size() - 1;i++){
            Node n = centroids[i];
            unsigned k = 0;
            unsigned dis = std::numeric_limits<unsigned>::max();
            for(unsigned j = i + 1;j < centroids.size();j++){
                unsigned tmp_dis = n.calculate_diff_kmp(centroids[j]);
                if(tmp_dis < dis){
                    dis = tmp_dis;
                    k = j;
                }
            }
            std::swap(centroids[i + 1], centroids[k]);
        }
        printf("---- [ 排序之后 ] ----:\n");
        for(unsigned i = 0;i < centroids.size();i++){
            printf("Center(%d): ", i);
            for(unsigned j = 0;j < centroids[i].feat.size();j++){
                printf("%u ", centroids[i].feat[j]);
            }
            printf("\n");
        }
        /* Kmeans算法执行 */
        unsigned int iter = 0;
        double cond = 0.08;
        double converge_rate = 0.1;
        while (iter < KMEANS_ITER) {
            // 清空聚类结果
            for (auto& cluster : clusters) 
                cluster.clear();
            // 将每个节点分配到最近的质心所在的聚类
            #pragma omp parallel for
            for (unsigned i = 0;i < nodes.size();i++) {
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { // 添加了关于尺寸的限制
                        min_diff = diff;
                        closestCentroid = j;
                    }
                }
                assert(closestCentroid != -1);
                nodes[i].cluster_id = closestCentroid;
            }
            #pragma omp parallel for
            for(unsigned i = 0; i < centroids.size();i++){
                for(unsigned nodeId = 0; nodeId < nodes.size(); nodeId++){
                    if(nodes[nodeId].cluster_id == i)
                        clusters[i].push_back(nodes[nodeId]);
                }
            }
            // 更新质心位置为聚类内节点的平均值
            unsigned num_not_converged = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < centroids.size(); ++i) {
                if (!clusters[i].empty()) {
                    unsigned dim = centroids[i].feat.size();
                    std::vector<unsigned> newCentroid(dim, 0);
                    for (size_t j = 0; j < dim; ++j){
                        for (const auto& node : clusters[i])
                            newCentroid[j] += node.feat[j];
                        double result = static_cast<double>(newCentroid[j]) / clusters[i].size();
                        newCentroid[j] =  static_cast<int>(std::round(result));
                    }
                    //unsigned dis = centroids[i].calculate_diff_vec(newCentroid);
                    unsigned dis = centroids[i].calculate_diff_vec(newCentroid);
                    if(dis >= 1){
                        #pragma omp atomic
                            num_not_converged++;
                    }
                    centroids[i].feat = newCentroid;
                }
            }
            iter++;
            printf("kmeans iter(%d)\n", iter);
            // 如果没有converge的分块占比少于10%，则结束算法
            if(num_not_converged < converge_rate * num_clusters){
                std::cout << "num not converged = " << num_not_converged << std::endl;
                break;
            }

        }
        for (unsigned i = 0;i < clusters.size();i++) 
            printf("cluster(%d) size = %d\n", i, clusters[i].size());
        // 直接得到new_id列表, 不做负载均衡
        int total_num = 0;
        int new_cnt = 0;
        for(unsigned cid = 0; cid < clusters.size();cid++){
            total_num += clusters[cid].size();
            for(unsigned nid = 0; nid < clusters[cid].size(); nid++){
                unsigned int old_id = clusters[cid][nid].id;
                new_id[old_id] = new_cnt;
                new_cnt++;
            }
        }
       printf("total_num = %d\n", total_num);
    }

    // [*] main hisorder algorithm (for push-pull mode)
    void HisOrder(){
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }
        std::vector<Node> centroids; // 质心数组
        unsigned num_clusters = this->graph->cluster_num;
        printf("分类数量 = %d\n", num_clusters);
        std::vector<std::vector<Node>> clusters(num_clusters); // 聚类数组

        /* KMeans++初始化 */
        int first_center_id = rand() % nodes.size();
        centroids.push_back(nodes[first_center_id]);
        for(unsigned i = 1; i < num_clusters;i++){
            if(i % 10  == 0)
                printf("Centroid %d init..\n", i);
            unsigned total_distance_sq = 0;
            std::vector<unsigned> distances(nodes.size(), std::numeric_limits<unsigned>::max());
            
            #pragma omp parallel for
            for(unsigned j = 0; j < nodes.size();j++){
                for(Node& centroid : centroids){
                    unsigned distance = centroid.calculate_diff_kmp(nodes[j]);
                    distances[j] = std::min(distances[j], distance);
                }
                unsigned sq = distances[j] * distances[j];
                #pragma omp atomic
                    total_distance_sq += sq;
            }
            double randValue = std::rand() / (RAND_MAX + 1.0) * total_distance_sq;
            for(unsigned j = 0; j < num_vertex;j++){
                randValue -= distances[j] * distances[j];
                if(randValue <= 0){
                    centroids.push_back(nodes[j]);
                    break;
                }
            }
        }
        printf("---- [ 排序之前 ] ----:\n");
        for(unsigned i = 0;i < centroids.size();i++){
            printf("Center(%d): ", i);
            for(unsigned j = 0;j < centroids[i].feat.size();j++){
                printf("%u ", centroids[i].feat[j]);
            }
            printf("\n");
        }
        for(unsigned i = 0;i < centroids.size() - 1;i++){
            Node n = centroids[i];
            unsigned k = 0;
            unsigned dis = std::numeric_limits<unsigned>::max();
            for(unsigned j = i + 1;j < centroids.size();j++){
                unsigned tmp_dis = n.calculate_diff_kmp(centroids[j]);
                if(tmp_dis < dis){
                    dis = tmp_dis;
                    k = j;
                }
            }
            std::swap(centroids[i + 1], centroids[k]);
        }
        printf("---- [ 排序之后 ] ----:\n");
        for(unsigned i = 0;i < centroids.size();i++){
            printf("Center(%d): ", i);
            for(unsigned j = 0;j < centroids[i].feat.size();j++){
                printf("%u ", centroids[i].feat[j]);
            }
            printf("\n");
        }
        /* KMeans算法执行 */
        unsigned int iter = 0;
        double cond = 0.08;
        double converge_rate = 0.1;
        while (iter < KMEANS_ITER) {
            // 清空聚类结果
            for (auto& cluster : clusters) 
                cluster.clear();
            // 将每个节点分配到最近的质心所在的聚类
            #pragma omp parallel for
            for (unsigned i = 0;i < nodes.size();i++) {
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { // 添加了关于尺寸的限制
                        min_diff = diff;
                        closestCentroid = j;
                    }
                }
                assert(closestCentroid != -1);
                nodes[i].cluster_id = closestCentroid;
            }
            #pragma omp parallel for
            for(unsigned i = 0; i < centroids.size();i++){
                for(unsigned nodeId = 0; nodeId < nodes.size(); nodeId++){
                    if(nodes[nodeId].cluster_id == i)
                        clusters[i].push_back(nodes[nodeId]);
                }
            }
            // 更新质心位置为聚类内节点的平均值
            unsigned num_not_converged = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < centroids.size(); ++i) {
                if (!clusters[i].empty()) {
                    unsigned dim = centroids[i].feat.size();
                    std::vector<unsigned> newCentroid(dim, 0);
                    for (size_t j = 0; j < dim; ++j){
                        for (const auto& node : clusters[i])
                            newCentroid[j] += node.feat[j];
                        double result = static_cast<double>(newCentroid[j]) / clusters[i].size();
                        newCentroid[j] =  static_cast<int>(std::round(result));
                    }
                    //unsigned dis = centroids[i].calculate_diff_vec(newCentroid);
                    unsigned dis = centroids[i].calculate_diff_vec(newCentroid);
                    if(dis >= 1){
                        #pragma omp atomic
                            num_not_converged++;
                    }
                    centroids[i].feat = newCentroid;
                }
            }
            iter++;
            printf("kmeans iter(%d)\n", iter);
            // 如果没有converge的分块占比少于10%，则结束算法
            if(num_not_converged < converge_rate * num_clusters){
                std::cout << "num not converged = " << num_not_converged << std::endl;
                break;
            }

        }
        for (unsigned i = 0;i < clusters.size();i++) 
            printf("cluster(%d) size = %d\n", i, clusters[i].size());
        /* 按照partition进行分块和负载均衡操作 */
        std::vector<Node> mapping;
        for(unsigned i = 0;i < clusters.size();i++){
            unsigned part_num = (clusters[i].size() - 1) / params::partition_size + 1;
            std::vector<std::vector<Node>> parts(part_num);
            // 构建large vertex和small vertex列表
            std::vector<Node> large_vertex;
            std::vector<Node> small_vertex;
            for(unsigned nid = 0; nid < clusters[i].size(); nid++){
                Node n = clusters[i][nid];
                if(graph->out_degree[n.id] > average_degree)
                    large_vertex.push_back(n);
                else
                    small_vertex.push_back(n);
            }
            unsigned p_id = 0;
            unsigned seg_size;
            /* 将large节点分段放进(8.18结束) */
            if(part_num > 1){
                for(unsigned lid = 0; lid < large_vertex.size();lid++){
                    if(parts[p_id].size() < params::partition_size){
                        parts[p_id++].push_back(large_vertex[lid]);
                        p_id = p_id % (part_num - 1);
                    }
                    else{
                        p_id = (p_id + 1) % (part_num - 1);
                        parts[part_num - 1].push_back(large_vertex[lid]);
                    }
                }
                // 将small节点成段放进
                seg_size = small_vertex.size() / part_num + 1;
                for(unsigned seg = 0; seg < part_num;seg++){
                    for(unsigned sid = seg * seg_size; sid < std::min((seg + 1) * seg_size, (unsigned)small_vertex.size());sid++){
                        if(parts[seg].size() < params::partition_size)
                            parts[seg].push_back(small_vertex[sid]);
                        else
                            parts[part_num - 1].push_back(small_vertex[sid]);
                    }
                }
            }
            else{
                seg_size = large_vertex.size();
                for(unsigned sid = 0;sid < seg_size;sid++){
                    parts[0].push_back(large_vertex[sid]);
                }
                seg_size = small_vertex.size();
                for(unsigned sid = 0;sid < seg_size;sid++){
                    parts[0].push_back(small_vertex[sid]);
                }
            }
            unsigned p = 0;
            for(p = 0; p < part_num;p++){
                for(unsigned k = 0; k < parts[p].size();k++)
                    mapping.push_back(parts[p][k]);
            }
        }
        printf("将所有节点分配完成后mapping数组的长度 = %ld\n", mapping.size());
        for(unsigned i = 0;i < mapping.size();i++){
            new_id[mapping[i].id] = i;
        }
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished\n");
        std::cout << "==========" << std::endl;

    }

    // [*] main hisorder algorithm (for push-only mode)
    void HO_SC() {
        /* 将入度为0的节点和普通节点区分开 */
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }
        std::vector<Node> centroids; // 质心数组
        unsigned num_clusters = this->graph->cluster_num;
        printf("分类数量 = %d\n", num_clusters);
        std::vector<std::vector<Node>> clusters(num_clusters); // 聚类数组

        /* KMeans++初始化 */
        int first_center_id = rand() % nodes.size();
        centroids.push_back(nodes[first_center_id]);
        for(unsigned i = 1; i < num_clusters;i++){
            if(i % 10  == 0)
                printf("Centroid %d init..\n", i);
            unsigned total_distance_sq = 0;
            std::vector<unsigned> distances(nodes.size(), std::numeric_limits<unsigned>::max());
            
            #pragma omp parallel for
            for(unsigned j = 0; j < nodes.size();j++){
                for(Node& centroid : centroids){
                    unsigned distance = centroid.calculate_diff_kmp(nodes[j]);
                    distances[j] = std::min(distances[j], distance);
                }
                unsigned sq = distances[j] * distances[j];
                #pragma omp atomic
                    total_distance_sq += sq;
            }
            double randValue = std::rand() / (RAND_MAX + 1.0) * total_distance_sq;
            for(unsigned j = 0; j < num_vertex;j++){
                randValue -= distances[j] * distances[j];
                if(randValue <= 0){
                    centroids.push_back(nodes[j]);
                    break;
                }
            }
        }
        printf("---- [ 排序之前 ] ----:\n");
        for(unsigned i = 0;i < centroids.size();i++){
            printf("Center(%d): ", i);
            for(unsigned j = 0;j < centroids[i].feat.size();j++){
                printf("%u ", centroids[i].feat[j]);
            }
            printf("\n");
        }
        for(unsigned i = 0;i < centroids.size() - 1;i++){
            Node n = centroids[i];
            unsigned k = 0;
            unsigned dis = std::numeric_limits<unsigned>::max();
            for(unsigned j = i + 1;j < centroids.size();j++){
                unsigned tmp_dis = n.calculate_diff_kmp(centroids[j]);
                if(tmp_dis < dis){
                    dis = tmp_dis;
                    k = j;
                }
            }
            std::swap(centroids[i + 1], centroids[k]);
        }
        printf("---- [ 排序之后 ] ----:\n");
        for(unsigned i = 0;i < centroids.size();i++){
            printf("Center(%d): ", i);
            for(unsigned j = 0;j < centroids[i].feat.size();j++){
                printf("%u ", centroids[i].feat[j]);
            }
            printf("\n");
        }
        /* KMeans算法执行 */
        unsigned int iter = 0;
        double cond = 0.08;
        double converge_rate = 0.1;
        while (iter < KMEANS_ITER) {
            // 清空聚类结果
            for (auto& cluster : clusters) 
                cluster.clear();
            // 将每个节点分配到最近的质心所在的聚类
            #pragma omp parallel for
            for (unsigned i = 0;i < nodes.size();i++) {
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    //unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { // 添加了关于尺寸的限制
                        min_diff = diff;
                        closestCentroid = j;
                    }
                }
                assert(closestCentroid != -1);
                nodes[i].cluster_id = closestCentroid;
            }
            #pragma omp parallel for
            for(unsigned i = 0; i < centroids.size();i++){
                for(unsigned nodeId = 0; nodeId < nodes.size(); nodeId++){
                    if(nodes[nodeId].cluster_id == i)
                        clusters[i].push_back(nodes[nodeId]);
                }
            }
            // 更新质心位置为聚类内节点的平均值
            unsigned num_not_converged = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < centroids.size(); ++i) {
                if (!clusters[i].empty()) {
                    unsigned dim = centroids[i].feat.size();
                    std::vector<unsigned> newCentroid(dim, 0);
                    for (size_t j = 0; j < dim; ++j){
                        for (const auto& node : clusters[i])
                            newCentroid[j] += node.feat[j];
                        double result = static_cast<double>(newCentroid[j]) / clusters[i].size();
                        newCentroid[j] =  static_cast<int>(std::round(result));
                    }
                    //unsigned dis = centroids[i].calculate_diff_vec(newCentroid);
                    unsigned dis = centroids[i].calculate_diff_vec(newCentroid);
                    if(dis >= 1){
                        #pragma omp atomic
                            num_not_converged++;
                        //#pragma omp critical
                            //std::cout << "diff = " << dis << std::endl;
                        
                    }
                    centroids[i].feat = newCentroid;
                }
            }
            iter++;
            printf("kmeans iter(%d)\n", iter);
            // 如果没有converge的分块占比少于10%，则结束算法
            if(num_not_converged < converge_rate * num_clusters){
                std::cout << "num not converged = " << num_not_converged << std::endl;
                break;
            }
        }
        for (unsigned i = 0;i < clusters.size();i++) 
            printf("cluster(%d) size = %d\n", i, clusters[i].size());
        /* 聚类结束, 开始获得partition聚类列表(考虑负载均衡) */
        std::vector<std::vector<Node>> parts(num_partitions);
        const auto average_degree = num_edges / num_vertex;
        int p_id = 0;
        int v_ptr = 0;
        for(unsigned cid = 0; cid < clusters.size();cid++){
            std::vector<Node> large_vertex;
            std::vector<Node> small_vertex;
            for(unsigned nid = 0; nid < clusters[cid].size();nid++){
                Node nd = clusters[cid][nid];
                if(graph->out_degree[nd.id] > average_degree)
                    large_vertex.push_back(nd);
                else
                    small_vertex.push_back(nd);
            }
            printf("large vertex(%u) = %ld\n", cid, large_vertex.size());
            // 将large节点全部放进
            for(unsigned lid = 0; lid < large_vertex.size();lid++){
                if(parts[p_id].size() < params::partition_size){
                    parts[p_id++].push_back(large_vertex[lid]);
                    p_id = p_id % (num_partitions - 1);
                }
                else{
                    p_id = (p_id + 1) % (num_partitions - 1);
                    parts[num_partitions - 1].push_back(large_vertex[lid]);
                }
            }
            // 将small节点成段放进
            unsigned seg_size = small_vertex.size() / (num_partitions - 1);
            for(unsigned seg = 0; seg < num_partitions - 1;seg++){
                for(unsigned sid = seg * seg_size; sid < (seg + 1) * seg_size;sid++){
                    if(parts[seg].size() < params::partition_size)
                        parts[seg].push_back(small_vertex[sid]);
                    else
                        parts[num_partitions - 1].push_back(small_vertex[sid]);
                }
            }
            for(unsigned sid = seg_size * (num_partitions - 1); sid < small_vertex.size();sid++){
                if(parts[p_id].size() < params::partition_size){
                    parts[p_id++].push_back(small_vertex[sid]);
                    p_id = p_id % (num_partitions - 1);
                }
                else{
                    p_id = (p_id + 1) % (num_partitions - 1);
                    parts[num_partitions - 1].push_back(small_vertex[sid]);
                }
            }
        }
        // 输出每个part的大小
        std::cout << "----------" << std::endl;
        for(unsigned i = 0;i < num_partitions;i++){
            printf("Part(%u): size = %d\n", i, parts[i].size());
        }
        unsigned int new_cnt = 0;
        for(unsigned pid = 0; pid < parts.size();pid++){
            for(unsigned nid = 0; nid < parts[pid].size(); nid++){
                unsigned int oid = parts[pid][nid].id;
                new_id[oid] = new_cnt;
                new_cnt++;
            }
        }
        //exit(1);
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
    }

    // 初始版本的horder
    void Hisorder_cc_sc() {
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }

        unsigned id = 0, pivot = 0, next_pivot = -1;
        unsigned cluster_id = 0;
        std::vector<std::vector<Node>> clusters; // 聚类数组
        
        while(id < num_vertex){
            std::vector<Node> cluster;
            clusters.push_back(cluster);
            Node pivot_node = nodes[pivot];
            next_pivot = -1;
            for(unsigned i = pivot; i < num_vertex;i++){
                if(nodes[i].cluster_id == -1 && pivot_node.calculate_diff(nodes[i]) == 0){
                    nodes[i].cluster_id = cluster_id;
                    clusters[cluster_id].push_back(nodes[i]);
                    id++;
                }
                else if(nodes[i].cluster_id == -1 && next_pivot == -1)
                    next_pivot = i; // 新建一个cluster, 并保存进clusters
            }
            cluster_id++;
            pivot = next_pivot;
        }
        printf("[Kmeans finish]: cluster num = %u\n", cluster_id);
        for (unsigned i = 0;i < clusters.size();i++) 
            printf("cluster(%d) size = %d\n", i, clusters[i].size());
        
        /* 聚类结束, 开始获得partition聚类列表(考虑负载均衡) */
        std::vector<std::vector<Node>> parts(num_partitions);
        const auto average_degree = num_edges / num_vertex;
        int p_id = 0;
        int v_ptr = 0;
        for(unsigned cid = 0; cid < clusters.size();cid++){
            std::vector<Node> large_vertex;
            std::vector<Node> small_vertex;
            for(unsigned nid = 0; nid < clusters[cid].size();nid++){
                Node nd = clusters[cid][nid];
                if(graph->out_degree[nd.id] > average_degree)
                    large_vertex.push_back(nd);
                else
                    small_vertex.push_back(nd);
            }
            printf("large vertex(%u) = %ld\n", cid, large_vertex.size());
            // 将large节点全部放进
            for(unsigned lid = 0; lid < large_vertex.size();lid++){
                if(parts[p_id].size() < params::partition_size){
                    parts[p_id++].push_back(large_vertex[lid]);
                    p_id = p_id % (num_partitions - 1);
                }
                else{
                    p_id = (p_id + 1) % (num_partitions - 1);
                    parts[num_partitions - 1].push_back(large_vertex[lid]);
                }
            }
            // 将small节点成段放进
            unsigned seg_size = small_vertex.size() / (num_partitions - 1);
            for(unsigned seg = 0; seg < num_partitions - 1;seg++){
                for(unsigned sid = seg * seg_size; sid < (seg + 1) * seg_size;sid++){
                    if(parts[seg].size() < params::partition_size)
                        parts[seg].push_back(small_vertex[sid]);
                    else
                        parts[num_partitions - 1].push_back(small_vertex[sid]);
                }
            }
            for(unsigned sid = seg_size * (num_partitions - 1); sid < small_vertex.size();sid++){
                if(parts[p_id].size() < params::partition_size){
                    parts[p_id++].push_back(small_vertex[sid]);
                    p_id = p_id % (num_partitions - 1);
                }
                else{
                    p_id = (p_id + 1) % (num_partitions - 1);
                    parts[num_partitions - 1].push_back(small_vertex[sid]);
                }
            }
        }
        // 输出每个part的大小
        std::cout << "----------" << std::endl;
        for(unsigned i = 0;i < num_partitions;i++){
            printf("Part(%u): size = %d\n", i, parts[i].size());
        }
        unsigned int new_cnt = 0;
        for(unsigned pid = 0; pid < parts.size();pid++){
            for(unsigned nid = 0; nid < parts[pid].size(); nid++){
                unsigned int oid = parts[pid][nid].id;
                new_id[oid] = new_cnt;
                new_cnt++;
            }
        }
        //exit(1);
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
    }

    // 初始版本的horder
    void HisOrder_cc() {
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }

        unsigned id = 0, pivot = 0, next_pivot = -1;
        unsigned cluster_id = 0;
        std::vector<std::vector<Node>> clusters; // 聚类数组
        
        while(id < num_vertex){
            std::vector<Node> cluster;
            clusters.push_back(cluster);
            Node pivot_node = nodes[pivot];
            next_pivot = -1;
            for(unsigned i = pivot; i < num_vertex;i++){
                if(nodes[i].cluster_id == -1 && pivot_node.calculate_diff(nodes[i]) == 0){
                    nodes[i].cluster_id = cluster_id;
                    clusters[cluster_id].push_back(nodes[i]);
                    id++;
                }
                else if(nodes[i].cluster_id == -1 && next_pivot == -1)
                    next_pivot = i; // 新建一个cluster, 并保存进clusters
            }
            cluster_id++;
            pivot = next_pivot;
        }
        printf("[Kmeans finish]: cluster num = %u\n", cluster_id);
        for (unsigned i = 0;i < clusters.size();i++) 
            printf("cluster(%d) size = %d\n", i, clusters[i].size());
        
        /* 按照partition进行分块和负载均衡操作 */
        std::vector<Node> mapping;
        for(unsigned i = 0;i < clusters.size();i++){
            unsigned part_num = (clusters[i].size() - 1) / params::partition_size + 1;
            std::vector<std::vector<Node>> parts(part_num);
            // 构建large vertex和small vertex列表
            std::vector<Node> large_vertex;
            std::vector<Node> small_vertex;
            for(unsigned nid = 0; nid < clusters[i].size(); nid++){
                Node n = clusters[i][nid];
                if(graph->out_degree[n.id] > average_degree)
                    large_vertex.push_back(n);
                else
                    small_vertex.push_back(n);
            }
            unsigned p_id = 0;
            unsigned seg_size;
            /* 将large节点分段放进(8.18结束) */
            if(part_num > 1){
                for(unsigned lid = 0; lid < large_vertex.size();lid++){
                    if(parts[p_id].size() < params::partition_size){
                        parts[p_id++].push_back(large_vertex[lid]);
                        p_id = p_id % (part_num - 1);
                    }
                    else{
                        p_id = (p_id + 1) % (part_num - 1);
                        parts[part_num - 1].push_back(large_vertex[lid]);
                    }
                }
                // 将small节点成段放进
                seg_size = small_vertex.size() / part_num + 1;
                for(unsigned seg = 0; seg < part_num;seg++){
                    for(unsigned sid = seg * seg_size; sid < std::min((seg + 1) * seg_size, (unsigned)small_vertex.size());sid++){
                        if(parts[seg].size() < params::partition_size)
                            parts[seg].push_back(small_vertex[sid]);
                        else
                            parts[part_num - 1].push_back(small_vertex[sid]);
                    }
                }
            }
            else{
                seg_size = large_vertex.size();
                for(unsigned sid = 0;sid < seg_size;sid++){
                    parts[0].push_back(large_vertex[sid]);
                }
                seg_size = small_vertex.size();
                for(unsigned sid = 0;sid < seg_size;sid++){
                    parts[0].push_back(small_vertex[sid]);
                }
            }
            unsigned p = 0;
            for(p = 0; p < part_num;p++){
                for(unsigned k = 0; k < parts[p].size();k++)
                    mapping.push_back(parts[p][k]);
            }
        }
        printf("将所有节点分配完成后mapping数组的长度 = %ld\n", mapping.size());
        for(unsigned i = 0;i < mapping.size();i++){
            new_id[mapping[i].id] = i;
        }
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished\n");
        std::cout << "==========" << std::endl;
    }

    // 没有负载均衡的cc版本
    void HisOrder_cc_noblc(){
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = -1;
        }
        unsigned id = 0, pivot = 0, next_pivot = -1;
        while(id < num_vertex){
            Node pivot_node = nodes[pivot];
            next_pivot = -1;
            for(unsigned i = pivot; i < num_vertex;i++){
                if(nodes[i].id == -1 && pivot_node.calculate_diff(nodes[i]) == 0){
                    nodes[i].id = id++;
                }
                else if(nodes[i].id == -1 && next_pivot == -1)
                    next_pivot = i;
            }
            pivot = next_pivot;
            next_pivot = -1;
        }
        printf("[Kmeans finish]: id = %u\n", id);

        // 得到new_id列表
        for(int i = 0; i < num_vertex;i++){
            new_id[i] = (unsigned)(nodes[i].id);
        }
    }

    // bfs序的history数据
    void HON_bfs() {
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].id = -1;
        }

        std::ifstream ifs(graph->in_feat); // 更改为实际的文件名
        if (!ifs.is_open()) {
            std::cout << "无法打开文件,程序退出" << std::endl;
            exit(1);
        }

        std::cout << "in feat = " << graph->in_feat << std::endl;

        unsigned vertexnum;
        ifs >> vertexnum;

        std::string line;
        unsigned nid = 0;
        while (std::getline(ifs, line)) {
            //std::cout << "This line = " << line << std::endl;
            std::vector<unsigned> row; // 用于存储每行的数据
            std::stringstream ss(line);
            std::string token;
            if(line == "----------")
                continue;
            while (std::getline(ss, token, ' ')) {
                unsigned value = (unsigned)(std::stoi(token));
                row.push_back(value);
            }
            std::sort(row.begin(), row.end());
            for(unsigned i = 0;i < row.size();i++){
                unsigned oid = row[i];
                if(nid == 0)
                    printf("start vertex = %u, new vertex = %u\n", oid, nid);
                nodes[oid].id = nid++;
            }
        }
        ifs.close();
    #ifdef DEBUG
        std::vector<unsigned> old_id(num_vertex);
    #endif
        // printf("start vertex = 606742, new vertex = %u\n", nodes[606742].id);
        for(int i = 0; i < num_vertex;i++){
            if(nodes[i].id == -1)
                nodes[i].id = nid++;
            new_id[i] = (unsigned)(nodes[i].id);
    #ifdef DEBUG 
            old_id[nodes[i].id] = i;
    #endif
        }
        std::cout << "final new id = " << nid << std::endl;
    #ifdef DEBUG
        std::ofstream outputFile("newid.txt");
        if (!outputFile.is_open()) {
            std::cout << "无法打开输出文件。" << std::endl;
            exit(1);
        }
        for (unsigned i = 0; i < num_vertex;i++) {
            outputFile << old_id[i] << std::endl;
        }
        outputFile.close();
    
    #endif
    }

    // bfs序的history数据, 并且实现负载均衡
    void HO_bfs() {
        /* 打开feat文件 */
        std::ifstream ifs(graph->in_feat); // 更改为实际的文件名
        if (!ifs.is_open()) {
            std::cout << "无法打开文件,程序退出" << std::endl;
            exit(1);
        }
        
        std::cout << graph->in_feat << std::endl;
        
        unsigned vertexnum;
        ifs >> vertexnum;
        
        /* 按照行读取特征文件 */
        std::string line;
        unsigned line_id = 0;
        std::vector<std::vector<unsigned>> cluster;
        std::vector<unsigned> visited(num_vertex, 0);
        while (std::getline(ifs, line)) {
            std::vector<unsigned> row; // 用于存储每行的数据
            std::stringstream ss(line);
            std::string token;
            if(line == "----------")
                continue;
            while (std::getline(ss, token, ' ')) {
                unsigned value = (unsigned)(std::stoi(token));
                visited[value] = 1;
                row.push_back(value);
            }
            if(row.size() > 0){
                std::sort(row.begin(), row.end());
                cluster.push_back(row);
                if(line_id == 0)
                    std::cout << "start vertex = " << row[0] << std::endl;
                line_id++;
            }
        }
        ifs.close();
        for(unsigned i = 0; i < cluster.size();i++){
            printf("cluster(%u): size = %d\n", i, cluster[i].size());
        }
        

        /* 聚类结束, 开始获得partition聚类列表(考虑负载均衡) */
        std::vector<std::vector<unsigned>> parts;
        const auto average_degree = num_edges / num_vertex;
        for(unsigned cid = 0; cid < cluster.size();cid++){
            printf("\n[cluster %d:]\n----------\n", cid);
            std::vector<unsigned> large_vertex;
            std::vector<unsigned> small_vertex;
            unsigned int old_id;
            for(unsigned nid = 0; nid < cluster[cid].size(); nid++){
                old_id = cluster[cid][nid];
                if(graph->out_degree[old_id] > average_degree)
                    large_vertex.push_back(old_id);
                else
                    small_vertex.push_back(old_id);
            }
            // 使用 shuffle 函数打乱 large vertex 中的元素            
            std::random_device rd;
            std::mt19937 rng(rd());
            std::shuffle(large_vertex.begin(), large_vertex.end(), rng);

            int num_partitions = (int)std::ceil((float)cluster[cid].size() / params::partition_size);
            printf("partition num = %d\n", num_partitions);
            int num_part_include_large = (num_partitions / params::num_threads) * params::num_threads;
            if(num_part_include_large == 0)
                num_part_include_large = num_partitions;
            int large_num = (int)std::ceil((float)large_vertex.size() / num_part_include_large);
            printf("num_part = %d, large_num = %d\n", num_part_include_large, large_num);
            // for every partition
            int sid = 0;
            for(unsigned part_id = 0; part_id < num_partitions;part_id++){
                std::vector<unsigned> part;
                long num_edge = 0;
                if(part_id < num_part_include_large){
                    for(int lid = part_id * large_num; lid < (part_id + 1) * large_num && lid < large_vertex.size();lid++){
                        part.push_back(large_vertex[lid]);
                        num_edge += graph->out_degree[large_vertex[lid]];
                    }
                    while(part.size() < params::partition_size && sid < small_vertex.size()){
                        num_edge += graph->out_degree[small_vertex[sid]];
                        part.push_back(small_vertex[sid++]);
                    }
                }
                else{
                    while(part.size() < params::partition_size && sid < small_vertex.size()){
                        num_edge += graph->out_degree[small_vertex[sid]];
                        part.push_back(small_vertex[sid++]);
                    }
                }
                parts.push_back(part);
                // printf("(%d:%ld:%ld) ", part_id, part.size(), num_edge);
            }
        }
        
        /* 另一种负载均衡方法: 随机shuffle 
        for(int i = 0;i < cluster.size();i++){
            std::random_device rd;
            std::mt19937 rng(rd());
            std::shuffle(cluster[i].begin(), cluster[i].end(), rng);
        }*/
        unsigned int new_cnt = 0;
        for(unsigned cid = 0; cid < cluster.size();cid++){
            for(unsigned nid = 0; nid < cluster[cid].size(); nid++){
                unsigned int oid = cluster[cid][nid];
                new_id[oid] = new_cnt;
                new_cnt++;
            }
        }
        // for vertices not visited yet
        for(unsigned i = 0;i < num_vertex;i++){
            if(visited[i] == 0){
                new_id[i] = new_cnt++;
            }
        }
        // for vertices not visited yet
        /*
        std::vector<unsigned> part;
        for(unsigned i = 0;i < num_vertex;i++){
            if(visited[i] == 0)
                part.push_back(i);
        }*/
        // 输出每个part的大小
        /*
        std::cout << "----------" << std::endl;
        for(unsigned i = 0;i < parts.size();i++){
            printf("Part(%u): size = %ld\n", i, parts[i].size());
        }
        */

        // 得到new_id列表
        /*
        for(unsigned pid = 0; pid < parts.size();pid++){
            for(unsigned nid = 0; nid < parts[pid].size(); nid++){
                unsigned int oid = parts[pid][nid];
                new_id[oid] = new_cnt;
                new_cnt++;
            }
        }
        */
        
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
        std::cout << "final new id = " << new_cnt << std::endl;
    }

    // bfs序的history数据, 并且根据线程的数量考虑负载均衡
    void HO_bfs_thread() {
        /* 通过feat文件构建cluster */
        unsigned vertexnum;
        std::ifstream ifs(graph->in_feat); // 更改为实际的文件名
        std::cout << graph->in_feat << std::endl;
        if (!ifs.is_open()) {
            std::cout << "无法打开文件,程序退出" << std::endl;
            exit(1);
        }
        ifs >> vertexnum;
        std::string line;
        unsigned line_id = 0;
        std::vector<std::vector<unsigned>> cluster;
        std::vector<unsigned> visited(num_vertex, 0);
        while (std::getline(ifs, line)) {
            std::vector<unsigned> row; // 用于存储每行的数据
            std::stringstream ss(line);
            std::string token;
            if(line == "----------")
                continue;
            while (std::getline(ss, token, ' ')) {
                unsigned value = (unsigned)(std::stoi(token));
                visited[value] = 1;
                row.push_back(value);
            }
            if(row.size() > 0){
                std::sort(row.begin(), row.end());
                cluster.push_back(row);
                if(line_id == 0)
                    std::cout << "start vertex = " << row[0] << std::endl;
                line_id++;
            }
        }
        ifs.close();
        std::cout << "----------" << std::endl;
        unsigned num_max_cluster = 0;
        for(unsigned i = 0; i < cluster.size();i++){
            if(cluster[i].size() > num_max_cluster)
                num_max_cluster = cluster[i].size();
            printf("cluster(%u): size = %d\n", i, cluster[i].size());
        }
        std::cout << "[!!] MAX cluster size = " << num_max_cluster << std::endl;
        #ifdef DEBUG
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据cluster信息, 构建chunk array列表 */
        unsigned num_thread = 48;
        unsigned col = cluster.size(); // 列数为cluster的数量
        unsigned row = (num_max_cluster / (params::partition_size * num_thread) + 1) * 48; // 行数为线程数量
        // 初始化chunk列表: 创建一个 n 行 m 列的二维 vector 数组，每个元素都是一个 vector<int>
        std::vector<std::vector<std::vector<int>>> chunk_array(row, std::vector<std::vector<int>>(col));
        const auto average_degree = num_edges / num_vertex;
        for(unsigned cid = 0; cid < cluster.size();cid++){
            unsigned row_id = 0;
            unsigned row_len = (cluster[cid].size() / (params::partition_size * 48) + 1) * 48;
            // 构建large vertex和small vertex列表
            std::vector<unsigned> large_vertex;
            std::vector<unsigned> small_vertex;
            for(unsigned nid = 0; nid < cluster[cid].size(); nid++){
                unsigned int v = cluster[cid][nid];
                if(graph->out_degree[v] > average_degree)
                    large_vertex.push_back(v);
                else
                    small_vertex.push_back(v);
            }
            // 将large节点全部放进
            for(unsigned lid = 0; lid < large_vertex.size();lid++){
                chunk_array[row_id++][cid].push_back(large_vertex[lid]);
                row_id = row_id % row_len;
            }
            // 将small节点成段放进
            unsigned seg_size = small_vertex.size() / row_len;
            for(unsigned seg_id = 0; seg_id < row_len;seg_id++){
                for(unsigned sid = seg_id * seg_size; sid < (seg_id + 1) * seg_size;sid++){
                    chunk_array[seg_id][cid].push_back(small_vertex[sid]);
                }
            }
            // 将剩余的small节点放进
            for(unsigned sid = seg_size * row_len; sid < small_vertex.size();sid++){
                chunk_array[row_id++][cid].push_back(small_vertex[sid]);
                row_id = row_id % row_len;
            }
        }
        #ifdef DEBUG
        std::cout << "-----------" << std::endl;
        std::cout << "展示chunk array:" << std::endl;
        for(unsigned i = 0;i < chunk_array.size();i++){
            printf("[ Row %u ]: ", i);
            for(unsigned j = 0;j < chunk_array[i].size();j++){
                std::cout << chunk_array[i][j].size() << "  ";
            }
            std::cout << std::endl;
        }
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据chunk_array, 构建part列表 */
        std::vector<std::vector<unsigned>> parts(num_partitions + 1);
        for(unsigned l = 0; l < row / num_thread;l++){
            for(unsigned i = 0; i < num_thread;i++){
                for(unsigned j = 0; j < col;j++){
                    auto chunk = chunk_array[num_thread * l + i][j];
                    if(chunk.size() > 0){
                        // 找到一个合适的partition: 尺寸合适并且id合适
                        unsigned p_id = i;
                        while(p_id < num_partitions && 
                                parts[p_id].size() + chunk.size() >= params::partition_size)
                            p_id += num_thread;
                        if(p_id < num_partitions){
                            parts[p_id].insert(parts[p_id].end(), chunk.begin(), chunk.end());
                        }
                        else{
                            parts[num_partitions].insert(parts[num_partitions].end(), chunk.begin(), chunk.end());
                        }
                    }
                    
                }
            }
        }
        // 将没有完整插入的分块按照节点的粒度进行插入
        for(unsigned i = 0;i < num_partitions - 1;i++){
            while(parts[i].size() < params::partition_size 
                    && parts[num_partitions - 1].size() > 0){
                unsigned v = parts[num_partitions - 1].back();
                parts[num_partitions - 1].pop_back();
                parts[i].push_back(v);
            }
        }
        // 将没有完整插入的分块按照节点的粒度进行插入
        for(unsigned i = 0;i < num_partitions;i++){
            while(parts[i].size() < params::partition_size 
                    && parts[num_partitions].size() > 0){
                unsigned v = parts[num_partitions].back();
                parts[num_partitions].pop_back();
                parts[i].push_back(v);
            }
        }
        // 将遍历中没有遍历到的节点插入到parts中去
        unsigned p_id = 0;
        for(unsigned i = 0;i < num_vertex;i++){
            if(visited[i] == 0){
                while(parts[p_id].size() >= params::partition_size)
                    p_id++;
                parts[p_id].push_back(i);
            }
        }
        // 输出结果
        std::cout << "----------" << std::endl;
        for(unsigned i = 0;i < num_partitions + 1;i++){
            printf("Part(%u): size = %d\n", i, parts[i].size());
        }
        #ifdef DEBUG
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据构造出来的partition列表, 得到new_id列表 */
        
        unsigned int new_cnt = 0;
        std::vector<unsigned> old_id(num_vertex, 0);
        for(unsigned pid = 0; pid < parts.size();pid++){
            for(unsigned nid = 0; nid < parts[pid].size(); nid++){
                unsigned int oid = parts[pid][nid];
                new_id[oid] = new_cnt;
                old_id[new_cnt] = oid;
                new_cnt++;
            }
        }
        //exit(1);
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
        std::cout << "final new id = " << new_cnt << std::endl;
    }

    // bfs序的history数据, 并且根据模式选择分块方法
    void HO_bfs_mode() {
        /* 通过feat文件构建cluster */
        unsigned vertexnum;
        std::ifstream ifs(graph->in_feat); // 更改为实际的文件名
        std::cout << graph->in_feat << std::endl;
        if (!ifs.is_open()) {
            std::cout << "无法打开文件,程序退出" << std::endl;
            exit(1);
        }
        ifs >> vertexnum;
        std::string line;
        unsigned line_id = 0;
        std::vector<std::vector<unsigned>> cluster;
        std::vector<unsigned> visited(num_vertex, 0);
        while (std::getline(ifs, line)) {
            std::vector<unsigned> row; // 用于存储每行的数据
            std::stringstream ss(line);
            std::string token;
            if(line == "----------")
                continue;
            while (std::getline(ss, token, ' ')) {
                unsigned value = (unsigned)(std::stoi(token));
                visited[value] = 1;
                row.push_back(value);
            }
            if(row.size() > 0){
                std::sort(row.begin(), row.end());
                cluster.push_back(row);
                if(line_id == 0)
                    std::cout << "start vertex = " << row[0] << std::endl;
                line_id++;
            }
        }
        ifs.close();
        std::cout << "----------" << std::endl;
        unsigned num_max_cluster = 0;
        for(unsigned i = 0; i < cluster.size();i++){
            if(cluster[i].size() > num_max_cluster)
                num_max_cluster = cluster[i].size();
            printf("cluster(%u): size = %d\n", i, cluster[i].size());
        }
        std::cout << "[!!] MAX cluster size = " << num_max_cluster << std::endl;
        #ifdef DEBUG
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据cluster信息, 构建chunk array列表 */
        unsigned num_thread = 48;
        unsigned col = cluster.size(); // 列数为cluster的数量
        unsigned row = (num_max_cluster / (params::partition_size * num_thread) + 1) * 48; // 行数为线程数量
        // 初始化chunk列表: 创建一个 n 行 m 列的二维 vector 数组，每个元素都是一个 vector<int>
        std::vector<std::vector<std::vector<int>>> chunk_array(row, std::vector<std::vector<int>>(col));
        const auto average_degree = num_edges / num_vertex;
        const auto average_num_edge = num_edges / num_partitions; // 每个part平均的边数
        std::vector<unsigned> pull_cluster_id;
        for(unsigned cid = 0; cid < cluster.size();cid++){
            unsigned row_id = 0;
            unsigned row_len = (cluster[cid].size() / (params::partition_size * 48) + 1) * 48;
            // 构建large vertex和small vertex列表
            std::vector<unsigned> large_vertex;
            std::vector<unsigned> small_vertex;
            unsigned cluster_num_edge = 0;
            for(unsigned nid = 0; nid < cluster[cid].size(); nid++){
                unsigned int v = cluster[cid][nid];
                cluster_num_edge += graph->out_degree[v];
                if(graph->out_degree[v] > average_degree)
                    large_vertex.push_back(v);
                else
                    small_vertex.push_back(v);
            }
            // 这个cluster使用Pull模式进行更新, 不再进行负载均衡
            if(cluster_num_edge > 0.5 * average_num_edge){
                pull_cluster_id.push_back(cid);
                continue;
            }
            // 将large节点全部放进
            for(unsigned lid = 0; lid < large_vertex.size();lid++){
                chunk_array[row_id++][cid].push_back(large_vertex[lid]);
                row_id = row_id % row_len;
            }
            // 将small节点成段放进
            unsigned seg_size = small_vertex.size() / row_len;
            for(unsigned seg_id = 0; seg_id < row_len;seg_id++){
                for(unsigned sid = seg_id * seg_size; sid < (seg_id + 1) * seg_size;sid++){
                    chunk_array[seg_id][cid].push_back(small_vertex[sid]);
                }
            }
            // 将剩余的small节点放进
            for(unsigned sid = seg_size * row_len; sid < small_vertex.size();sid++){
                chunk_array[row_id++][cid].push_back(small_vertex[sid]);
                row_id = row_id % row_len;
            }
        }
        #ifdef DEBUG
        std::cout << "-----------" << std::endl;
        std::cout << "展示chunk array:" << std::endl;
        for(unsigned i = 0;i < chunk_array.size();i++){
            printf("[ Row %u ]: ", i);
            for(unsigned j = 0;j < chunk_array[i].size();j++){
                std::cout << chunk_array[i][j].size() << "  ";
            }
            std::cout << std::endl;
        }
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        std::vector<std::vector<unsigned>> parts(num_partitions + 1);
        // 先插入pull模式的节点
        unsigned pid = 0;
        for(unsigned i = 0;i < pull_cluster_id.size();i++){
            unsigned cluster_id = pull_cluster_id[i];
            for(unsigned j = 0;j < cluster[cluster_id].size();j++){
                unsigned vertex = cluster[cluster_id][j];
                if(parts[pid].size() >= params::partition_size){
                    pid++;
                }
                parts[pid].push_back(vertex);
            }
        }
        #ifdef DEBUG
        std::cout << "----------" << std::endl;
        std::cout << "Pull模式节点更新完成" << std::endl;
        for(unsigned i = 0;i <= pid;i++){
            printf("Parts(%u): size = %u\n", i, parts[i].size());
        }
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据chunk_array, 构建part列表 */
        for(unsigned l = 0; l < row / num_thread;l++){
            for(unsigned i = 0; i < num_thread;i++){
                for(unsigned j = 0; j < col;j++){
                    auto chunk = chunk_array[num_thread * l + i][j];
                    if(chunk.size() > 0){
                        // 找到一个合适的partition: 尺寸合适并且id合适
                        unsigned p_id = i;
                        while(p_id < num_partitions && 
                                parts[p_id].size() + chunk.size() >= params::partition_size)
                            p_id += num_thread;
                        if(p_id < num_partitions){
                            parts[p_id].insert(parts[p_id].end(), chunk.begin(), chunk.end());
                        }
                        else{
                            parts[num_partitions].insert(parts[num_partitions].end(), chunk.begin(), chunk.end());
                        }
                    }
                    
                }
            }
        }
        // 将没有完整插入的分块按照节点的粒度进行插入
        for(unsigned i = 0;i < num_partitions - 1;i++){
            while(parts[i].size() < params::partition_size 
                    && parts[num_partitions - 1].size() > 0){
                unsigned v = parts[num_partitions - 1].back();
                parts[num_partitions - 1].pop_back();
                parts[i].push_back(v);
            }
        }
        // 将没有完整插入的分块按照节点的粒度进行插入
        for(unsigned i = 0;i < num_partitions;i++){
            while(parts[i].size() < params::partition_size 
                    && parts[num_partitions].size() > 0){
                unsigned v = parts[num_partitions].back();
                parts[num_partitions].pop_back();
                parts[i].push_back(v);
            }
        }
        // 将遍历中没有遍历到的节点插入到parts中去
        unsigned p_id = 0;
        for(unsigned i = 0;i < num_vertex;i++){
            if(visited[i] == 0){
                while(parts[p_id].size() >= params::partition_size)
                    p_id++;
                parts[p_id].push_back(i);
            }
        }
        // 输出结果
        std::cout << "----------" << std::endl;
        for(unsigned i = 0;i < num_partitions + 1;i++){
            printf("Part(%u): size = %d\n", i, parts[i].size());
        }
        #ifdef DEBUG
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据构造出来的partition列表, 得到new_id列表 */
        unsigned int new_cnt = 0;
        std::vector<unsigned> old_id(num_vertex, 0);
        for(unsigned pid = 0; pid < parts.size();pid++){
            for(unsigned nid = 0; nid < parts[pid].size(); nid++){
                unsigned int oid = parts[pid][nid];
                new_id[oid] = new_cnt;
                old_id[new_cnt] = oid;
                new_cnt++;
            }
        }
        //exit(1);
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
        std::cout << "final new id = " << new_cnt << std::endl;
    }

    // 考虑极度冷节点的history-based grouping
    void HO_cold() {
        unsigned num_vertex_norm = num_vertex - graph->num_vertex_deg_0 - graph->num_vertex_deg_1;
        std::vector<Node> deg0_node(graph->num_vertex_deg_0);
        std::vector<Node> deg1_node(graph->num_vertex_deg_1);
        std::vector<Node> norm_node(num_vertex_norm);
        //std::vector<Node> norm_node(num_vertex);
        // 初始化三个数组
        unsigned deg0_cnt, deg1_cnt, norm_cnt;
        deg0_cnt = deg1_cnt = norm_cnt = 0;
        for(unsigned i = 0;i < num_vertex;i++){
            if(graph->in_degree[i] == 0){
                deg0_node[deg0_cnt].id = i;
                deg0_node[deg0_cnt].cluster_id = -1;
                deg0_cnt++;
            }
            else if(graph->in_degree[i] > 0){
                deg1_node[deg1_cnt].id = i;
                deg1_node[deg1_cnt].cluster_id = -1;
                deg1_node[deg1_cnt].src_id = graph->in_degree[i] - 1;
                deg1_cnt++;
            }
            else{
                norm_node[norm_cnt].id = i;
                norm_node[norm_cnt].cluster_id = -1;
                norm_node[norm_cnt].feat = graph->attr[i];
                norm_cnt++;
            }
        }
        std::cout << "[Value Check]: " << std::endl;
        std::cout << "\t#v(deg = 0) = " << deg0_cnt << std::endl;
        std::cout << "\t#v(deg = 1) = " << deg1_cnt << std::endl;
        std::cout << "\t#v(deg > 1) = " << norm_cnt << std::endl;
        std::cout << "-----------" << std::endl;
        // Kmeans算法初始化
        double cond = 0.1;
        double converge_rate = 0.1;
        std::vector<Node> centroids; // 质心数组
        std::vector<std::vector<Node>> clusters(num_partitions); // 聚类数组

        // 随机初始化
        for(int i = 0;i < num_partitions;i++){
            int rand_node = rand() % norm_node.size();
            centroids.push_back(norm_node[rand_node]);
        }

        unsigned int iter = 0;
        while (iter < KMEANS_ITER) {
            // 清空聚类结果
            for (auto& cluster : clusters) 
                cluster.clear();
            
            // 将每个节点分配到最近的质心所在的聚类
            #pragma omp parallel for
            for (unsigned i = 0;i < norm_node.size();i++) {
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    //unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    unsigned diff = norm_node[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { // 添加了关于尺寸的限制
                        min_diff = diff;
                        closestCentroid = j;
                    }
                }
                assert(closestCentroid != -1);
                norm_node[i].cluster_id = closestCentroid;
            }

            #pragma omp parallel for
            for(unsigned i = 0; i < centroids.size();i++){
                for(unsigned nodeId = 0; nodeId < norm_node.size(); nodeId++){
                    if(norm_node[nodeId].cluster_id == i)
                        clusters[i].push_back(norm_node[nodeId]);
                }
            }

            // 更新质心位置为聚类内节点的平均值
            unsigned num_not_converged = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < centroids.size(); ++i) {
                if (!clusters[i].empty()) {
                    unsigned dim = centroids[i].feat.size();
                    std::vector<unsigned> newCentroid(dim, 0);
                    for (size_t j = 0; j < dim; ++j){
                        for (const auto& node : clusters[i])
                            newCentroid[j] += node.feat[j];
                        double result = static_cast<double>(newCentroid[j]) / clusters[i].size();
                        newCentroid[j] =  static_cast<int>(std::round(result));
                    }
                    if(centroids[i].calculate_diff_vec(newCentroid) > cond * newCentroid.size()){
                        #pragma omp atomic
                            num_not_converged++;
                    }
                    centroids[i].feat = newCentroid;
                }
            }
            printf("iter(%d): not converged = (%d / %d), converge rate = %lf\n", 
                    iter, num_not_converged, num_partitions, 
                    (double)(num_partitions - num_not_converged) / num_partitions);
            iter++;
            // 如果没有converge的分块占比少于10%，则结束算法
            if(num_not_converged < converge_rate * num_partitions){
                printf("iter(End): not converged = (%d / %d)\n", num_not_converged, num_partitions);
                break;
            }
        }
        std::vector<unsigned> cluster_size(clusters.size());
        for (unsigned i = 0;i < clusters.size();i++) {
            cluster_size[i] = clusters[i].size();
            // printf("cluster(%d) size = %d\n", i, clusters[i].size());
        }
        std::sort(cluster_size.begin(), cluster_size.end(), [](int a, int b) { return a > b; });
        std::cout << "----------" << std::endl;
        int num = 0;
        for(unsigned i = 0;i < cluster_size.size();i++){
            num += cluster_size[i];
            printf("Cluster(%d) -> %u\n", i, cluster_size[i]);
        }
        // 聚类结束
        // 得到new_id列表
        std::vector<int> n2c(num_vertex, 0);
        int new_cnt = 0;
        for(unsigned cid = 0; cid < clusters.size();cid++){
            for(unsigned nid = 0; nid < clusters[cid].size(); nid++){
                unsigned int old_id = clusters[cid][nid].id;
                n2c[old_id] = cid + 1;
                new_id[old_id] = new_cnt;
                new_cnt++;
            }
        }
        // 对冷节点进行操作
        unsigned num_src_noc = 0;
        for(unsigned i = 0;i < deg1_node.size();i++){
            if(n2c[deg1_node[i].src_id] != 0){
                deg1_node[i].cluster_id = n2c[deg1_node[i].src_id] - 1;
            }
            else{
                deg1_node[i].cluster_id = std::numeric_limits<int>::max();
                num_src_noc++;
            }
        }
        std::sort(deg1_node.begin(), deg1_node.end(), [](Node a, Node b) { return a.cluster_id < b.cluster_id; });
        // 为度数等于1的节点赋予id
        for(unsigned i = 0;i < deg1_node.size();i++){
            new_id[deg1_node[i].id] = new_cnt++;
        }
        // 为度数等于0的节点赋予新id
        for(unsigned i = 0;i < deg0_cnt;i++){
            new_id[deg0_node[i].id] = new_cnt++;
        }
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
        // exit(1);
    }

    // 考虑极度冷节点的history-based grouping
    void HO_cold_blc() {
        unsigned num_vertex_norm = num_vertex - graph->num_vertex_deg_0 - graph->num_vertex_deg_1;
        std::vector<Node> deg0_node(graph->num_vertex_deg_0);
        std::vector<Node> deg1_node(graph->num_vertex_deg_1);
        std::vector<Node> norm_node(num_vertex_norm);
        //std::vector<Node> norm_node(num_vertex);
        // 初始化三个数组
        unsigned deg0_cnt, deg1_cnt, norm_cnt;
        deg0_cnt = deg1_cnt = norm_cnt = 0;
        for(unsigned i = 0;i < num_vertex;i++){
            if(graph->in_degree[i] == 0){
                deg0_node[deg0_cnt].id = i;
                deg0_node[deg0_cnt].cluster_id = -1;
                deg0_cnt++;
            }
            else if(graph->in_degree[i] > 0){
                deg1_node[deg1_cnt].id = i;
                deg1_node[deg1_cnt].cluster_id = -1;
                deg1_node[deg1_cnt].src_id = graph->in_degree[i] - 1;
                deg1_cnt++;
            }
            else{
                norm_node[norm_cnt].id = i;
                norm_node[norm_cnt].cluster_id = -1;
                norm_node[norm_cnt].feat = graph->attr[i];
                norm_cnt++;
            }
        }
        std::cout << "[Value Check]: " << std::endl;
        std::cout << "\t#v(deg = 0) = " << deg0_cnt << std::endl;
        std::cout << "\t#v(deg = 1) = " << deg1_cnt << std::endl;
        std::cout << "\t#v(deg > 1) = " << norm_cnt << std::endl;
        std::cout << "-----------" << std::endl;
        // Kmeans算法初始化
        double cond = 0.1;
        double converge_rate = 0.1;
        std::vector<Node> centroids; // 质心数组
        std::vector<std::vector<Node>> clusters(num_partitions); // 聚类数组

        // 随机初始化
        for(int i = 0;i < num_partitions;i++){
            int rand_node = rand() % norm_node.size();
            centroids.push_back(norm_node[rand_node]);
        }

        unsigned int iter = 0;

        while (iter < KMEANS_ITER) {
            // 清空聚类结果
            for (auto& cluster : clusters) 
                cluster.clear();
            
            // 将每个节点分配到最近的质心所在的聚类
            #pragma omp parallel for
            for (unsigned i = 0;i < norm_node.size();i++) {
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    //unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    unsigned diff = norm_node[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { // 添加了关于尺寸的限制
                        min_diff = diff;
                        closestCentroid = j;
                    }
                }
                assert(closestCentroid != -1);
                norm_node[i].cluster_id = closestCentroid;
            }

            #pragma omp parallel for
            for(unsigned i = 0; i < centroids.size();i++){
                for(unsigned nodeId = 0; nodeId < norm_node.size(); nodeId++){
                    if(norm_node[nodeId].cluster_id == i)
                        clusters[i].push_back(norm_node[nodeId]);
                }
            }

            // 更新质心位置为聚类内节点的平均值
            unsigned num_not_converged = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < centroids.size(); ++i) {
                if (!clusters[i].empty()) {
                    unsigned dim = centroids[i].feat.size();
                    std::vector<unsigned> newCentroid(dim, 0);
                    for (size_t j = 0; j < dim; ++j){
                        for (const auto& node : clusters[i])
                            newCentroid[j] += node.feat[j];
                        double result = static_cast<double>(newCentroid[j]) / clusters[i].size();
                        newCentroid[j] =  static_cast<int>(std::round(result));
                    }
                    if(centroids[i].calculate_diff_vec(newCentroid) > cond * newCentroid.size()){
                        #pragma omp atomic
                            num_not_converged++;
                    }
                    centroids[i].feat = newCentroid;
                }
            }
            printf("iter(%d): not converged = (%d / %d), converge rate = %lf\n", 
                    iter, num_not_converged, num_partitions, 
                    (double)(num_partitions - num_not_converged) / num_partitions);
            iter++;
            // 如果没有converge的分块占比少于10%，则结束算法
            if(num_not_converged < converge_rate * num_partitions){
                printf("iter(End): not converged = (%d / %d)\n", num_not_converged, num_partitions);
                break;
            }
        }
        std::vector<unsigned> cluster_size(clusters.size());
        for (unsigned i = 0;i < clusters.size();i++) {
            cluster_size[i] = clusters[i].size();
            // printf("cluster(%d) size = %d\n", i, clusters[i].size());
        }
        std::sort(cluster_size.begin(), cluster_size.end(), [](int a, int b) { return a > b; });
        std::cout << "----------" << std::endl;
        int num = 0;
        for(unsigned i = 0;i < cluster_size.size();i++){
            num += cluster_size[i];
            printf("Cluster(%d) -> %u\n", i, cluster_size[i]);
        }

        /* 聚类结束, 开始获得partition聚类列表(考虑负载均衡) */
        std::vector<std::vector<unsigned>> parts(num_partitions);
        const auto average_degree = num_edges / num_vertex;
        // 得到new_id列表
        
        int p_id = 0;
        for(unsigned cid = 0; cid < clusters.size();cid++){
            std::vector<unsigned> large_vertex;
            std::vector<unsigned> small_vertex;
            unsigned int old_id;
            for(unsigned nid = 0; nid < clusters[cid].size(); nid++){
                old_id = clusters[cid][nid].id;
                if(graph->out_degree[old_id] > average_degree)
                    large_vertex.push_back(old_id);
                else
                    small_vertex.push_back(old_id);
            }
            // 将large节点全部放进
            for(unsigned lid = 0; lid < large_vertex.size();lid++){
                parts[p_id++].push_back(large_vertex[lid]);
                p_id = p_id % num_partitions;
            }
            // 将small节点成段放进
            unsigned seg_size = small_vertex.size() / num_partitions;
            for(unsigned seg = 0; seg < num_partitions;seg++){
                for(unsigned sid = seg * seg_size; sid < (seg + 1) * seg_size;sid++){
                    parts[seg].push_back(small_vertex[sid]);
                }
            }
            for(unsigned sid = seg_size * num_partitions; sid < small_vertex.size();sid++){
                parts[p_id++].push_back(small_vertex[sid]);
                p_id = p_id % num_partitions;
            }
        }
        // 分块结束
        std::vector<int> n2c(num_vertex, 0);
        for(unsigned cid = 0; cid < clusters.size();cid++){
            for(unsigned nid = 0; nid < clusters[cid].size(); nid++){
                unsigned int old_id = clusters[cid][nid].id;
                n2c[old_id] = cid + 1;
            }
        }
        // 对冷节点进行操作
        unsigned num_src_noc = 0;
        for(unsigned i = 0;i < deg1_node.size();i++){
            if(n2c[deg1_node[i].src_id] != 0){
                deg1_node[i].cluster_id = n2c[deg1_node[i].src_id] - 1;
            }
            else{
                deg1_node[i].cluster_id = std::numeric_limits<int>::max();
                num_src_noc++;
            }
        }
        std::sort(deg1_node.begin(), deg1_node.end(), [](Node a, Node b) { return a.cluster_id < b.cluster_id; });
        // 对于度数为1的节点
        unsigned i = 0;
        for(unsigned part_id = 0; part_id < num_partitions;part_id++){
            while(i < deg1_node.size() && 
                    parts[part_id].size() < num_vertex / num_partitions + 1){
                parts[part_id].push_back(deg1_node[i].id);
                i++;
            }
        }
        // 对于度数为0的节点
        i = 0;
        for(unsigned part_id = 0; part_id < num_partitions;part_id++){
            while(i < deg0_node.size() && 
                    parts[part_id].size() < num_vertex / num_partitions + 1){
                parts[part_id].push_back(deg0_node[i].id);
                i++;
            }
        }
        
        // 得到new_id列表
        int new_cnt = 0;
        for(unsigned pid = 0; pid < parts.size();pid++){
            for(unsigned nid = 0; nid < parts[pid].size(); pid++){
                unsigned int old_id = parts[pid][nid];
                new_id[old_id] = new_cnt;
                new_cnt++;
            }
        }
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
        // exit(1);
    }

    // kmeans序的history数据, 根据模式选择分块方法
    void HO_mode() {
        /* 将入度为0的节点和普通节点区分开 */
        std::vector<Node> nodes; // 节点数组
        std::vector<unsigned> deg_zero;
        for(unsigned i = 0;i < num_vertex;i++){
            if(graph->in_degree[i] == 0){
                deg_zero.push_back(i);
            }
            else{
                Node nd;
                nd.feat = graph->attr[i];
                nd.id = i;
                nd.cluster_id = -1;
                nodes.push_back(nd);
            }
        }
        printf("num deg_zero = %ld\n", deg_zero.size());
        //exit(1);
        /* 根据特征向量, 对入度不为0的节点进行聚类 */
        double cond = 0.08;
        double converge_rate = 0.01;
        std::vector<Node> centroids; // 质心数组
        // unsigned num_cluster = this->graph->cluster_num;
        unsigned num_cluster = params::num_partitions;
        std::vector<std::vector<unsigned>> cluster(num_cluster); // 聚类数组
        for(int i = 0;i < num_cluster;i++){
            int rand_node = rand() % nodes.size(); // 随机初始化
            centroids.push_back(nodes[rand_node]);
        }
        unsigned int iter = 0;
        while (iter < KMEANS_ITER) {
            for (auto& c : cluster) 
                c.clear();// 清空聚类结果
            #pragma omp parallel for
            for (unsigned i = 0;i < nodes.size();i++) { // 将每个节点分配到最近的质心所在的聚类
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { // 添加了关于尺寸的限制
                        min_diff = diff;
                        closestCentroid = j;
                    }
                }
                assert(closestCentroid != -1);
                nodes[i].cluster_id = closestCentroid;
            }
            #pragma omp parallel for
            for(unsigned i = 0; i < centroids.size();i++){
                for(unsigned nodeId = 0; nodeId < nodes.size(); nodeId++){
                    if(nodes[nodeId].cluster_id == i)
                        cluster[i].push_back(nodeId); // 对节点聚类
                }
            }
            unsigned num_not_converged = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < centroids.size(); ++i) {// 更新质心位置为聚类内节点的平均值
                if (!cluster[i].empty()) {
                    unsigned dim = centroids[i].feat.size();
                    std::vector<unsigned> newCentroid(dim, 0);
                    for (size_t j = 0; j < dim; ++j){
                        for (const auto& nodeId : cluster[i])
                            newCentroid[j] += nodes[nodeId].feat[j];
                        double result = static_cast<double>(newCentroid[j]) / cluster[i].size();
                        newCentroid[j] = static_cast<int>(std::round(result));
                    }
                    if(centroids[i].calculate_diff_vec(newCentroid) > cond * newCentroid.size()){
                        #pragma omp atomic
                            num_not_converged++;
                    }
                    centroids[i].feat = newCentroid;
                }
            }
            printf("iter(%d): not converged = (%d / %d), converge rate = %lf\n", 
                    iter, num_not_converged, num_cluster, 
                    (double)(num_cluster - num_not_converged) / num_cluster);
            iter++;
            if(num_not_converged < converge_rate * num_cluster){
                printf("iter(End): not converged = (%d / %d)\n", num_not_converged, num_cluster);
                break;// 如果没有converge的分块占比少于10%，则结束算法
            }
        }
        printf("----------\n");
        unsigned num_max_cluster = 0;
        for(unsigned i = 0; i < cluster.size();i++){ // 输出聚类结果
            if(cluster[i].size() > num_max_cluster)
                num_max_cluster = cluster[i].size();
            for(unsigned j = 0;j < cluster[i].size();j++){
                cluster[i][j] = nodes[cluster[i][j]].id;
            }
            printf("center(%u): ", i);
            for(unsigned k = 0;k < centroids[i].feat.size();k++)
                printf("%u ", centroids[i].feat[k]);
            printf("\n");
            printf("cluster(%u): size = %ld\n", i, cluster[i].size());
        }
        std::cout << "[!!] MAX cluster size = " << num_max_cluster << std::endl;
        #ifdef DEBUG
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据cluster信息, 构建chunk array列表 */
        unsigned num_thread = 48;
        unsigned col = cluster.size(); // 列数为cluster的数量
        unsigned row = (num_max_cluster / (params::partition_size * num_thread) + 1) * 48; // 行数为线程数量
        // 初始化chunk列表: 创建一个 n 行 m 列的二维 vector 数组，每个元素都是一个 vector<int>
        std::vector<std::vector<std::vector<int>>> chunk_array(row, std::vector<std::vector<int>>(col));
        const auto average_degree = num_edges / num_vertex;
        const auto average_num_edge = num_edges / num_partitions; // 每个part平均的边数
        std::vector<unsigned> pull_cluster_id;
        for(unsigned cid = 0; cid < cluster.size();cid++){
            unsigned row_id = 0;
            unsigned row_len = (cluster[cid].size() / (params::partition_size * 48) + 1) * 48;
            // 构建large vertex和small vertex列表
            std::vector<unsigned> large_vertex;
            std::vector<unsigned> small_vertex;
            unsigned cluster_num_edge = 0;
            for(unsigned nid = 0; nid < cluster[cid].size(); nid++){
                unsigned int v = cluster[cid][nid];
                cluster_num_edge += graph->out_degree[v];
                if(graph->out_degree[v] > average_degree)
                    large_vertex.push_back(v);
                else
                    small_vertex.push_back(v);
            }
            // 这个cluster使用Pull模式进行更新, 不再进行负载均衡
            if(cluster_num_edge > 0.5 * average_num_edge){
                pull_cluster_id.push_back(cid);
                continue;
            }
            // 将large节点全部放进
            for(unsigned lid = 0; lid < large_vertex.size();lid++){
                chunk_array[row_id++][cid].push_back(large_vertex[lid]);
                row_id = row_id % row_len;
            }
            // 将small节点成段放进
            unsigned seg_size = small_vertex.size() / row_len;
            for(unsigned seg_id = 0; seg_id < row_len;seg_id++){
                for(unsigned sid = seg_id * seg_size; sid < (seg_id + 1) * seg_size;sid++){
                    chunk_array[seg_id][cid].push_back(small_vertex[sid]);
                }
            }
            // 将剩余的small节点放进
            for(unsigned sid = seg_size * row_len; sid < small_vertex.size();sid++){
                chunk_array[row_id++][cid].push_back(small_vertex[sid]);
                row_id = row_id % row_len;
            }
        }
        #ifdef DEBUG
        std::cout << "-----------" << std::endl;
        std::cout << "展示chunk array:" << std::endl;
        for(unsigned i = 0;i < chunk_array.size();i++){
            printf("[ Row %u ]: ", i);
            for(unsigned j = 0;j < chunk_array[i].size();j++){
                std::cout << chunk_array[i][j].size() << "  ";
            }
            std::cout << std::endl;
        }
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        printf("----------\n");
        printf("所有pull-based cluster的ID: ");
        for(unsigned i = 0;i < pull_cluster_id.size();i++)
            printf("%u ", pull_cluster_id[i]);
        printf("\n");
        std::vector<std::vector<unsigned>> parts(num_partitions + 1);
        /* 将pull模式的节点全部插入 */
        unsigned pid = 0;
        for(unsigned i = 0;i < pull_cluster_id.size();i++){
            unsigned cluster_id = pull_cluster_id[i];
            for(unsigned j = 0;j < cluster[cluster_id].size();j++){
                unsigned vertex = cluster[cluster_id][j];
                if(parts[pid].size() >= params::partition_size){
                    pid++;
                }
                parts[pid].push_back(vertex);
            }
        }
        #ifdef DEBUG
        std::cout << "----------" << std::endl;
        std::cout << "Pull模式节点更新完成" << std::endl;
        for(unsigned i = 0;i <= pid;i++){
            printf("Parts(%u): size = %u\n", i, parts[i].size());
        }
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据chunk_array, 构建part列表 */
        for(unsigned l = 0; l < row / num_thread;l++){
            for(unsigned i = 0; i < num_thread;i++){
                for(unsigned j = 0; j < col;j++){
                    auto chunk = chunk_array[num_thread * l + i][j];
                    if(chunk.size() > 0){
                        // 找到一个合适的partition: 尺寸合适并且id合适
                        unsigned p_id = i;
                        while(p_id < num_partitions && 
                                parts[p_id].size() + chunk.size() >= params::partition_size)
                            p_id += num_thread;
                        if(p_id < num_partitions){
                            parts[p_id].insert(parts[p_id].end(), chunk.begin(), chunk.end());
                        }
                        else{
                            parts[num_partitions].insert(parts[num_partitions].end(), chunk.begin(), chunk.end());
                        }
                    }
                    
                }
            }
        }
        for(unsigned i = 0;i < num_partitions - 1;i++){
            while(parts[i].size() < params::partition_size // 将没有完整插入的分块按照节点的粒度进行插入
                    && parts[num_partitions - 1].size() > 0){
                unsigned v = parts[num_partitions - 1].back();
                parts[num_partitions - 1].pop_back();
                parts[i].push_back(v);
            }
        }
        for(unsigned i = 0;i < num_partitions;i++){
            while(parts[i].size() < params::partition_size 
                    && parts[num_partitions].size() > 0){
                unsigned v = parts[num_partitions].back();
                parts[num_partitions].pop_back();
                parts[i].push_back(v);// 将没有完整插入的分块按照节点的粒度进行插入
            }
        }
        // 将遍历中没有遍历到的节点插入到parts中去
        unsigned p_id = 0;
        for(auto v: deg_zero){
            while(parts[p_id].size() >= params::partition_size)
                p_id++;
            parts[p_id].push_back(v);
        }
        // 输出结果
        std::cout << "----------" << std::endl;
        for(unsigned i = 0;i < num_partitions + 1;i++){
            printf("Part(%u): size = %ld\n", i, parts[i].size());
        }
        #ifdef DEBUG
        std::cout << "按  Enter  键继续" << std::endl;
        std::cin.get();
        #endif
        /* 根据构造出来的partition列表, 得到new_id列表 */
        unsigned int new_cnt = 0;
        std::vector<unsigned> old_id(num_vertex, 0);
        for(unsigned pid = 0; pid < parts.size();pid++){
            for(unsigned nid = 0; nid < parts[pid].size(); nid++){
                unsigned int oid = parts[pid][nid];
                new_id[oid] = new_cnt;
                old_id[new_cnt] = oid;
                new_cnt++;
            }
        }
        //exit(1);
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
    }

    // 初始版本的horder
    void HisOrder_it() {
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }
        std::vector<Node> centroids; // 质心数组
        unsigned num_clusters = this->graph->cluster_num;
        printf("分类数量 = %d\n", num_clusters);
        std::vector<std::vector<Node>> clusters(num_clusters); // 聚类数组

        /* KMeans++初始化 */
        int first_center_id = rand() % nodes.size();
        centroids.push_back(nodes[first_center_id]);
        for(unsigned i = 1; i < num_clusters;i++){
            if(i % 10  == 0)
                printf("Centroid %d init..\n", i);
            unsigned total_distance_sq = 0;
            std::vector<unsigned> distances(num_vertex, std::numeric_limits<unsigned>::max());
            
            #pragma omp parallel for
            for(unsigned j = 0; j < nodes.size();j++){
                for(Node& centroid : centroids){
                    unsigned distance = centroid.calculate_diff_kmp(nodes[j]);
                    distances[j] = std::min(distances[j], distance);
                }
                unsigned sq = distances[j] * distances[j];
                #pragma omp atomic
                    total_distance_sq += sq;
            }
            double randValue = std::rand() / (RAND_MAX + 1.0) * total_distance_sq;
            for(unsigned j = 0; j < num_vertex;j++){
                randValue -= distances[j] * distances[j];
                if(randValue <= 0){
                    centroids.push_back(nodes[j]);
                    break;
                }
            }
        }
        printf("---- [ 排序之前 ] ----:\n");
        for(unsigned i = 0;i < centroids.size();i++){
            printf("Center(%d): ", i);
            for(unsigned j = 0;j < centroids[i].feat.size();j++){
                printf("%u ", centroids[i].feat[j]);
            }
            printf("\n");
        }
        for(unsigned i = 0;i < centroids.size() - 1;i++){
            Node n = centroids[i];
            unsigned k = 0;
            unsigned dis = std::numeric_limits<unsigned>::max();
            for(unsigned j = i + 1;j < centroids.size();j++){
                unsigned tmp_dis = n.calculate_diff_kmp(centroids[j]);
                if(tmp_dis < dis){
                    dis = tmp_dis;
                    k = j;
                }
            }
            std::swap(centroids[i + 1], centroids[k]);
        }
        printf("---- [ 排序之后 ] ----:\n");
        for(unsigned i = 0;i < centroids.size();i++){
            printf("Center(%d): ", i);
            for(unsigned j = 0;j < centroids[i].feat.size();j++){
                printf("%u ", centroids[i].feat[j]);
            }
            printf("\n");
        }
        /* Kmeans算法执行 */
        unsigned int iter = 0;
        double cond = 0.08;
        double converge_rate = 0.1;
        while (iter < KMEANS_ITER) {
            // 清空聚类结果
            for (auto& cluster : clusters) 
                cluster.clear();
            // 将每个节点分配到最近的质心所在的聚类
            #pragma omp parallel for
            for (unsigned i = 0;i < nodes.size();i++) {
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { // 添加了关于尺寸的限制
                        min_diff = diff;
                        closestCentroid = j;
                    }
                }
                assert(closestCentroid != -1);
                nodes[i].cluster_id = closestCentroid;
            }
            #pragma omp parallel for
            for(unsigned i = 0; i < centroids.size();i++){
                for(unsigned nodeId = 0; nodeId < nodes.size(); nodeId++){
                    if(nodes[nodeId].cluster_id == i)
                        clusters[i].push_back(nodes[nodeId]);
                }
            }
            // 更新质心位置为聚类内节点的平均值
            unsigned num_not_converged = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < centroids.size(); ++i) {
                if (!clusters[i].empty()) {
                    unsigned dim = centroids[i].feat.size();
                    std::vector<unsigned> newCentroid(dim, 0);
                    for (size_t j = 0; j < dim; ++j){
                        for (const auto& node : clusters[i])
                            newCentroid[j] += node.feat[j];
                        double result = static_cast<double>(newCentroid[j]) / clusters[i].size();
                        newCentroid[j] =  static_cast<int>(std::round(result));
                    }
                    //unsigned dis = centroids[i].calculate_diff_vec(newCentroid);
                    unsigned dis = centroids[i].calculate_diff_vec(newCentroid);
                    if(dis >= 1){
                        #pragma omp atomic
                            num_not_converged++;
                    }
                    centroids[i].feat = newCentroid;
                }
            }
            iter++;
            printf("kmeans iter(%d)\n", iter);
            // 如果没有converge的分块占比少于10%，则结束算法
            if(num_not_converged < converge_rate * num_clusters){
                std::cout << "num not converged = " << num_not_converged << std::endl;
                break;
            }

        }
        for (unsigned i = 0;i < clusters.size();i++) 
            printf("cluster(%d) size = %d\n", i, clusters[i].size());
        /* 聚类结束, 开始获得partition聚类列表(考虑负载均衡) */
        std::vector<std::vector<unsigned>> parts;
        const auto average_degree = num_edges / num_vertex;
        for(unsigned cid = 0; cid < clusters.size();cid++){
            printf("\n[cluster %d:]\n----------\n", cid);
            std::vector<unsigned> large_vertex;
            std::vector<unsigned> small_vertex;
            unsigned int old_id;
            for(unsigned nid = 0; nid < clusters[cid].size(); nid++){
                old_id = clusters[cid][nid].id;
                if(graph->in_degree[old_id] > average_degree)
                    large_vertex.push_back(old_id);
                else
                    small_vertex.push_back(old_id);
            }
            // 使用 shuffle 函数打乱 large vertex 中的元素            
            std::random_device rd;
            std::mt19937 rng(rd());
            std::shuffle(large_vertex.begin(), large_vertex.end(), rng);

            int num_partitions = (int)std::ceil((float)clusters[cid].size() / params::partition_size);
            printf("partition num = %d\n", num_partitions);
            int num_part_include_large = (num_partitions / params::num_threads) * params::num_threads;
            if(num_part_include_large == 0)
                num_part_include_large = num_partitions;
            int large_num = (int)std::ceil((float)large_vertex.size() / num_part_include_large);
            printf("num_part = %d, large_num = %d\n", num_part_include_large, large_num);
            // for every partition
            int sid = 0;
            for(unsigned part_id = 0; part_id < num_partitions;part_id++){
                std::vector<unsigned> part;
                long num_edge = 0;
                if(part_id < num_part_include_large){
                    for(int lid = part_id * large_num; lid < (part_id + 1) * large_num && lid < large_vertex.size();lid++){
                        part.push_back(large_vertex[lid]);
                        num_edge += graph->in_degree[large_vertex[lid]];
                    }
                    while(part.size() < params::partition_size && sid < small_vertex.size()){
                        num_edge += graph->in_degree[small_vertex[sid]];
                        part.push_back(small_vertex[sid++]);
                    }
                }
                else{
                    while(part.size() < params::partition_size && sid < small_vertex.size()){
                        num_edge += graph->in_degree[small_vertex[sid]];
                        part.push_back(small_vertex[sid++]);
                    }
                }
                parts.push_back(part);
                printf("(%d:%ld:%ld) ", part_id, part.size(), num_edge);
            }
        }
        // 输出每个part的大小
        std::cout << "----------" << std::endl;
        for(unsigned i = 0;i < parts.size();i++){
            printf("Part(%u): size = %ld\n", i, parts[i].size());
        }

        // 得到new_id列表
        int new_cnt = 0;
        for(unsigned pid = 0; pid < parts.size();pid++){
            for(unsigned nid = 0; nid < parts[pid].size(); nid++){
                unsigned int oid = parts[pid][nid];
                new_id[oid] = new_cnt;
                new_cnt++;
            }
        }
        
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished, new count = %d\n", new_cnt);
        std::cout << "==========" << std::endl;
        std::cout << "final new id = " << new_cnt << std::endl;
    }

    /* 根据重排获得的new_id列表, 获得新的出度列表/csr数据(row列表和col列表) */
    void getNewGraph(Algo algo) {
        cpu_timer timer;
        float time = 0.0;
        unsigned max_threads = omp_get_max_threads();

        std::vector<unsigned> new_degree(num_vertex, 0);
        timer.start();
        //Assign the outdegree to new id
        #pragma omp parallel for schedule(static) num_threads(max_threads)
        for(unsigned i = 0; i < num_vertex; i++){
            new_degree[new_id[i]] = graph->out_degree[i];
        }
        float tm = timer.elapsed().wall/(1e9);

        // Build new row_index array
        std::vector<unsigned> new_row(num_vertex + 1, 0);
        __gnu_parallel::partial_sum(new_degree.begin(), new_degree.end(), new_row.begin() + 1);
        std::vector<unsigned> new_col(num_edges, 0);
        tm = timer.elapsed().wall/(1e9); 

        #ifdef WEIGHTED
         std::vector<unsigned> new_wei(num_edges, 0);
        #endif
        //Build new col_index array
        #pragma omp parallel for schedule(dynamic) num_threads(max_threads)
        for(unsigned i = 0; i < num_vertex; i++) {
            unsigned count = 0;
            for(unsigned j = graph->row_index[i]; j < graph->row_index[i + 1]; j++) {
                new_col[new_row[new_id[i]] + count] = new_id[graph->col_index[j]];
                #ifdef WEIGHTED
                new_wei[new_row[new_id[i]] + count] = graph->edge_weight[j];
                #endif
                count++;
            }
            std::sort(new_col.begin() + new_row[new_id[i]], new_col.begin() + new_row[new_id[i]] + count);
        }
        tm = timer.elapsed().wall/(1e9); 

        this->graph->out_degree.swap(new_degree);
        this->graph->row_index.swap(new_row);
        this->graph->col_index.swap(new_col);

        #ifdef WEIGHTED
        new_graph.edge_weight.swap(new_wei);
        #endif
    }

    /* 重排算法: 目标是按照给定的算法获得new_id列表 */
    void reorder(Algo algo, int vertex){
        switch(algo) {
            case Algo::original: 
                std::cout << "[!!] original order is maintained" << '\n';
                break;
            case Algo::hisorder_wo_blc:
                std::cout << "[!!] reordering method: [ history-based order ]" << '\n';
                HisOrder_wo_blc();
                break;
            case Algo::hisorder:
                std::cout << "[!!] reordering method: [ history-based order + balance ]" << '\n';
                HisOrder();
                break;
            case Algo::HO_SC:
                std::cout << "[!!] Reordering method: [kmeans + balance + SC]" << std::endl;
                HO_SC();
                break;
            case Algo::Hisorder_cc_sc:
                std::cout << "[!!] Reordering method: [cc order + sc mode]" << std::endl;
                Hisorder_cc_sc();
                break;
            case Algo::HisOrder_cc:
                std::cout << "[!!] Reordering method: [cc order + mix mode]" << std::endl;
                HisOrder_cc();
                break;
            case Algo::HisOrder_cc_noblc:
                std::cout << "[!!] Reordering method: [cc order + no balance]" << std::endl;
                HisOrder_cc_noblc();
                break;
            case Algo::HON_bfs:
                std::cout << "[!!] Reordering method: [bfs order]" << std::endl;
                HON_bfs();
                break;
            case Algo::HO_bfs:
                std::cout << "[!!] Reordering method: [bfs order + balance]" << std::endl;
                HO_bfs();
                break;
            case Algo::bfs_blc:
                std::cout << "[!!] Reordering method: [bfs order + balance + consider thread + mode choice]" << std::endl;
                HO_bfs_thread();
                break;
            /*
            case Algo::HO_bfs_thread:
                std::cout << "[!!] Reordering method: [bfs order + balance + consider thread]" << std::endl;
                HO_bfs_thread();
                break;
            
            case Algo::HON_cold:
                std::cout << "[!!] Reordering method: [kmeans + cold]" << std::endl;
                HO_cold();
                break;
            case Algo::HO_cold:
                std::cout << "[!!] Reordering method: [kmeans + cold + balance]" << std::endl;
                HO_cold_blc();
                break;
            case Algo::HO_mode:
                std::cout << "[!!] Reordering method: [kmeans + balance + consider thread + mode choice]" << std::endl;
                HO_mode();
                break;
            */
            
            /*
            case Algo::hisorder_it:
                std::cout << "[!!] Reordering method: [kmeans + balance for graphit]" << std::endl;
                HisOrder_it();
                break;
            */
            default:
                std::cout << "choose a correct algorithm!" << '\n';
        }
        std::cout << "[NEW]:  start vertex id: " << new_id[vertex] << '\n';
    }
};