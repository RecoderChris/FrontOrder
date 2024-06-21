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

class Node {
public:
    int id;
    int cluster_id;
    unsigned src_id;
    std::vector<unsigned> feat;

    unsigned calculate_diff(Node node) {
        unsigned diff = 0;
        for (size_t i = 0; i < feat.size(); i++){
            int m = (int)feat[i] - (int)node.feat[i];
            diff += std::abs(m);
        }
        return diff;
    }

    unsigned calculate_diff_kmp(Node node){
        unsigned diff = 0;
        for (size_t i = 0; i < feat.size(); i++){
            long int m = (long int)feat[i] - (long int)node.feat[i];
            diff += std::abs(m);
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

    unsigned calculate_diff_vec(std::vector<unsigned> feat_vec){
        unsigned diff = 0;
        for(size_t i = 0; i < feat.size();i++){
            int m = (int)feat[i] - (int)feat_vec[i];
            diff += std::abs(m);
    
        }
        return diff;
    }

};

bool compareByValue(const Node& a, const Node& b) {
    return a.cluster_id < b.cluster_id;
}

class Orderer {
    typedef std::pair<unsigned, unsigned> degree_id_pair; 
    unsigned num_vertex;
    unsigned num_edges; 
    unsigned num_partitions; 
    unsigned num_levels; 
    unsigned average_degree; 
    std::vector<unsigned> levels; 

    Graph* graph;
public:
    std::vector<unsigned> new_id;
    std::vector<unsigned> new_dst_id;

    Orderer(Graph* g) {
        this->graph = g;
        num_vertex = graph->num_vertex;
        num_edges = graph->num_edges;
        num_partitions = params::num_partitions;
        new_id = std::vector<unsigned>(num_vertex, 0);
        new_dst_id = std::vector<unsigned>(num_vertex, 0);
        average_degree = num_edges/num_vertex;

        num_levels = (unsigned) log2((float) num_partitions) + 2;

        for(int i = 0; i < num_levels; i++)
            levels.push_back(average_degree * pow(2, i - 1));

        levels.back() = UINT_MAX;
    }

    void HisOrder_wo_blc() {
        std::vector<Node> nodes(num_vertex); 
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }
        std::vector<Node> centroids; 
        unsigned num_clusters = this->graph->cluster_num;
        printf("分类数量 = %d\n", num_clusters);
        std::vector<std::vector<Node>> clusters(num_clusters); 

        /* KMeans++ */
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
        unsigned int iter = 0;
        double cond = 0.08;
        double converge_rate = 0.1;
        while (iter < KMEANS_ITER) {
            for (auto& cluster : clusters) 
                cluster.clear();
            #pragma omp parallel for
            for (unsigned i = 0;i < nodes.size();i++) {
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { 
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
            if(num_not_converged < converge_rate * num_clusters){
                std::cout << "num not converged = " << num_not_converged << std::endl;
                break;
            }

        }
        for (unsigned i = 0;i < clusters.size();i++) 
            printf("cluster(%d) size = %d\n", i, clusters[i].size());
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
    void HisOrder_PCPM(){
        std::vector<Node> nodes(num_vertex); 
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }
        std::vector<Node> centroids; 
        unsigned num_clusters = this->graph->cluster_num;
        printf("分类数量 = %d\n", num_clusters);
        std::vector<std::vector<Node>> clusters(num_clusters); 

        /* KMeans++初始�? */
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
        unsigned int iter = 0;
        double cond = 0.08;
        double converge_rate = 0.1;
        while (iter < KMEANS_ITER) {
            for (auto& cluster : clusters) 
                cluster.clear();
            #pragma omp parallel for
            for (unsigned i = 0;i < nodes.size();i++) {
                unsigned min_diff = std::numeric_limits<unsigned>::max();
                int closestCentroid = -1;
                for (int j = 0; j < centroids.size(); j++) {
                    unsigned diff = nodes[i].calculate_diff(centroids[j]);
                    if (diff < min_diff) { 
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
            if(num_not_converged < converge_rate * num_clusters){
                std::cout << "num not converged = " << num_not_converged << std::endl;
                break;
            }

        }
        for (unsigned i = 0;i < clusters.size();i++) 
            printf("cluster(%d) size = %d\n", i, clusters[i].size());
        std::vector<Node> mapping;
        for(unsigned i = 0;i < clusters.size();i++){
            unsigned part_num = (clusters[i].size() - 1) / params::partition_size + 1;
            std::vector<std::vector<Node>> parts(part_num);
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
        printf("将所有节点分配完成后mapping数组的长�? = %ld\n", mapping.size());
        for(unsigned i = 0;i < mapping.size();i++){
            new_id[mapping[i].id] = i;
        }
        std::cout << "----------" << std::endl;
        printf("[!!]Reorder finished\n");
        std::cout << "==========" << std::endl;

    }

    // [*] main hisorder algorithm (for push-only mode)
    void HisOrder() {
        std::vector<Node> nodes(num_vertex); // 节点数组
        #pragma omp parallel for
        for(unsigned i = 0;i < num_vertex;i++){
            nodes[i].feat = graph->attr[i];
            nodes[i].id = i;
            nodes[i].cluster_id = -1;
        }
        std::vector<Node> centroids; 
        unsigned num_clusters = this->graph->cluster_num;
        printf("分类数量 = %d\n", num_clusters);
        std::vector<std::vector<Node>> clusters(num_clusters); // 聚类数组

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
            // 更新质心位置为聚类内节点的平均�?
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
            // 如果没有converge的分块占比少�?10%，则结束算法
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
        // 输出每个part的大�?
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

    /* Random */
    void fastRandom() {        
        std::vector<unsigned> index(num_vertex);
        #pragma omp parallel for
        for(unsigned i = 0; i < num_vertex; i++)
            index[i] = i; 

        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(index), std::end(index), rng);

        #pragma omp parallel for
        for(unsigned i = 0; i < num_vertex; i++)
            new_id[index[i]] = i;  

    }

    /* Sorting by Degree */
    void fastSort() {
        std::vector<unsigned>& out_degree = graph->out_degree;
        std::vector<degree_id_pair> degree_vs_id(num_vertex);

        #pragma omp parallel for
        for(unsigned i = 0; i < num_vertex; i++)
            degree_vs_id[i] = std::make_pair(out_degree[i], i);

        __gnu_parallel::sort(degree_vs_id.begin(), degree_vs_id.end(), std::greater<degree_id_pair>());
        
        std::cout << degree_vs_id.front().first << " " << degree_vs_id.front().second << '\n';

        #pragma omp parallel for
        for(unsigned i = 0; i < num_vertex; i++)
            new_id[degree_vs_id[i].second] = i;
    }

    /* Degree-based Group */
    void fastDBG(unsigned num_levels) {

        levels.clear();
        for(int i = 0; i < num_levels; i++)
            levels.push_back(average_degree * pow(2, i - 1));
        levels.back() = UINT_MAX;

        const auto& out_degree = graph->out_degree;

        std::vector<unsigned>segment[params::num_threads][num_levels];
        
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned i = 0; i < num_vertex; i++) {
            for(unsigned j = 0; j < num_levels; j++) {
                if(out_degree[i] <= levels[j]) {
                    segment[omp_get_thread_num()][j].push_back(i);
                    break;
                }
            }
        }
        unsigned tmp = 0;
        unsigned seg_offset[params::num_threads][num_levels];
        for(int j = num_levels - 1; j >= 0; j--)
            for(unsigned t = 0; t < params::num_threads; t++) {
                seg_offset[t][j] = tmp;
                tmp += segment[t][j].size();
            }
        
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned t = 0; t < params::num_threads; t++)
            for(int j = num_levels - 1; j >= 0; j--) {
                unsigned offset = seg_offset[t][j];
                const std::vector<unsigned>& curr_seg = segment[t][j];
                for(auto id: curr_seg)
                    new_id[id] = offset++;
            }
    }

    /* Frequency-based Cluster */
    void fastFBC() {
        new_id.clear();
        new_id = std::move(std::vector<unsigned>(num_vertex, UINT_MAX));

        Vector2d<degree_id_pair> local_degree_vs_id(params::num_threads);
        unsigned slice = num_vertex / params::num_threads;
        std::vector<unsigned> start(params::num_threads); 
        std::vector<unsigned> end(params::num_threads);
        std::vector<unsigned> hub_count(params::num_threads);
        std::vector<unsigned> non_hub_count(params::num_threads);
        std::vector<unsigned> index(params::num_threads);
        std::vector<unsigned>& out_degree = graph->out_degree;
        unsigned sum_hub_count = 0;

        for(int t = 0; t < params::num_threads; t++) {
            start[t] = t * slice;
            end[t] = (t + 1) * slice;
            hub_count[t] = 0;
        }
        end[params::num_threads - 1] = num_vertex;

        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned t = 0; t < params::num_threads; t++) 
            for(unsigned i = start[t]; i < end[t]; i++)
                if(out_degree[i] > average_degree) 
                    local_degree_vs_id[t].push_back(std::make_pair(out_degree[i], i));

        std::vector<unsigned> hub_offset(params::num_threads + 1, 0);
        for(int t = 0; t < params::num_threads; t++) {
            hub_count[t] = local_degree_vs_id[t].size();
            sum_hub_count += hub_count[t];
            non_hub_count[t] = end[t] - start[t] - hub_count[t];
            hub_offset[t + 1] =  hub_offset[t] + hub_count[t];
        }
        index[0] = sum_hub_count;
        for(int t = 1; t < params::num_threads; t++) 
            index[t] = index[t - 1] + non_hub_count[t - 1];

        std::vector<degree_id_pair> degree_vs_id(sum_hub_count);

        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(int i = 0; i < params::num_threads; i++) {
            for(unsigned j = 0; j < hub_count[i]; j++)
                degree_vs_id[hub_offset[i]++] = local_degree_vs_id[i][j];
            std::vector<degree_id_pair>().swap(local_degree_vs_id[i]);
        }
        __gnu_parallel::sort(degree_vs_id.begin(), degree_vs_id.end(),
            std::greater<degree_id_pair>());

        #pragma omp parallel for
        for(unsigned i = 0; i < sum_hub_count; i++)
            new_id[degree_vs_id[i].second] = i;
        
        std::vector<degree_id_pair>().swap(degree_vs_id);

        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(int t = 0; t < params::num_threads; t++)
            for(unsigned i = start[t]; i < end[t]; i++)
                if(new_id[i] == UINT_MAX) 
                    new_id[i] = index[t]++;

    }

    /* HubCluster */
    void fastHC() {
        levels.clear();
        num_levels = 2;
        levels.push_back(average_degree);
        levels.push_back(UINT_MAX);

        const auto& out_degree = graph->out_degree;

        std::vector<unsigned>segment[params::num_threads][num_levels];
        
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned i = 0; i < num_vertex; i++) {
            for(unsigned j = 0; j < num_levels; j++) {
                if(out_degree[i] <= levels[j]) {
                    segment[omp_get_thread_num()][j].push_back(i);
                    break;
                }
            }
        }
        unsigned tmp = 0;
        unsigned seg_offset[params::num_threads][num_levels];
        for(int j = num_levels - 1; j >= 0; j--)
            for(unsigned t = 0; t < params::num_threads; t++) {
                seg_offset[t][j] = tmp;
                tmp += segment[t][j].size();
            }
        
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned t = 0; t < params::num_threads; t++)
            for(int j = num_levels - 1; j >= 0; j--) {
                unsigned offset = seg_offset[t][j];
                const std::vector<unsigned>& curr_seg = segment[t][j];
                for(auto id: curr_seg)
                    new_id[id] = offset++;
            }
    }

    /* Cache-aware Reorder */
    void fastCorder() {
        unsigned max_threads = omp_get_max_threads();

        Vector2d<unsigned> large_segment(max_threads);
        Vector2d<unsigned> small_segment(max_threads);

        const auto average_degree = num_edges/num_vertex;
        
        // classifying hot/cold vertices
        // the static scheduler ensures the relative order 
        // (static, 1024) disrupts the relative order but achieves very good performance
        #pragma omp parallel for schedule(static) num_threads(max_threads) 
        for(unsigned i = 0; i < num_vertex; i++)
            if(graph->out_degree[i] > average_degree) 
                large_segment[omp_get_thread_num()].push_back(i);
            else
                small_segment[omp_get_thread_num()].push_back(i);

        std::vector<unsigned> large_offset(max_threads + 1, 0);
        std::vector<unsigned> small_offset(max_threads + 1, 0);

        large_offset[1] = large_segment[0].size();
        small_offset[1] = small_segment[0].size(); 
        for(unsigned i = 0; i < max_threads ; i++) {
            large_offset[i+1] = large_offset[i] + large_segment[i].size();
            small_offset[i+1] = small_offset[i] + small_segment[i].size();
        }

       unsigned total_large = large_offset[max_threads];
       unsigned total_small = small_offset[max_threads]; 
     
        unsigned num_large_per_seg = ceil((float) total_large  / num_partitions);
        unsigned num_small_per_seg = params::partition_size - num_large_per_seg;

        unsigned last_cls = num_partitions - 1;

      //  constructing partitions based on the classified hot/cold vertices
        #pragma omp parallel for schedule(static) num_threads(max_threads)
        for(unsigned i = 0; i < num_partitions; i++) {
            unsigned index = i * params::partition_size;
            unsigned num_large =  (i != num_partitions - 1) ? (i + 1) * num_large_per_seg: total_large;
            unsigned large_start_t = 0;
            unsigned large_end_t = 0;
            unsigned large_start_v = 0;
            unsigned large_end_v = 0;
            unsigned large_per_seg = (i != num_partitions - 1) ? num_large_per_seg: total_large - i * num_large_per_seg;

            unsigned num_small =  (i != num_partitions - 1) ? (i + 1) * num_small_per_seg: total_small;
            unsigned small_start_t = 0;
            unsigned small_end_t = 0;
            unsigned small_start_v = 0;
            unsigned small_end_v = 0;
            unsigned small_per_seg = (i != num_partitions - 1) ? num_small_per_seg: total_small - i * num_small_per_seg;
            //HOT find the starting segment and starting vertex
            for(unsigned t = 0; t < max_threads; t++) {
                if(large_offset[t+1] > num_large - large_per_seg) {
                    large_start_t = t;
                    large_start_v = num_large - large_per_seg - large_offset[t];
                    break;
                }
            }
            //HOT find the ending segment and ending vertex
            for(unsigned t = large_start_t; t < max_threads; t++) {
                if(large_offset[t+1] >= num_large) {
                    large_end_t = t;
                    large_end_v =  num_large - large_offset[t] - 1;
                    break;
                }
            }

            //COLD find the starting segment and starting vertex
            for(unsigned t = 0; t < max_threads; t++) {
                if(small_offset[t+1] > num_small - small_per_seg) {
                    small_start_t = t;
                    small_start_v = num_small - small_per_seg - small_offset[t];
                    break;
                }
            }
            //COLD find the ending segment and ending vertex
           for(unsigned t = small_start_t; t < max_threads; t++) {
                if(small_offset[t+1] >= num_small) {
                    small_end_t = t;
                    small_end_v =  num_small - small_offset[t] - 1;
                    break;
                }
            }
            // HOT move the vertices form hot segment(s) to a partition
            if(large_start_t == large_end_t) {
                for(unsigned j = large_start_v; j <= large_end_v; j++) {
                    new_id[large_segment[large_start_t][j]] = index++;
                }
            } else {
                for(unsigned t = large_start_t; t < large_end_t; t++) {
                    if(t!=large_start_t)
                        large_start_v = 0;
                    for(unsigned j = large_start_v; j < large_segment[t].size(); j++) {
                        new_id[large_segment[t][j]] = index++;
                    }
                }
                for(unsigned j = 0; j <= large_end_v; j++) {
                    new_id[large_segment[large_end_t][j]] = index++;
                }
            }
            // COLD move the vertices form cold segment(s) to a partition
            if(small_start_t == small_end_t) {
                for(unsigned j = small_start_v; j <= small_end_v; j++) {
                    new_id[small_segment[small_start_t][j]] = index++;
                }
            } else {
                for(unsigned t = small_start_t; t < small_end_t; t++) {
                    if(t!=small_start_t)
                        small_start_v = 0;
                    for(unsigned j = small_start_v; j < small_segment[t].size(); j++) {
                        new_id[small_segment[t][j]] = index++;
                    }
                }
                for(unsigned j = 0; j <= small_end_v; j++) {
                    new_id[small_segment[small_end_t][j]] = index++;
                }
            }
        }
    }

    // mapping graph reorder
    void mapReorder(std::string mapping_file) {
        std::ifstream ifs(mapping_file.c_str(), std::ifstream::in);
        if (!ifs.good()) {
            std::cout << "File " << mapping_file << " does not exist!" << std::endl;
            exit(-1);
        }
        unsigned long int num_vertex_1, num_edges_1;
        ifs >> num_vertex_1;
        ifs >> num_edges_1;
        std::cout << " num_vertex: " << num_vertex_1 << " num_edges: "  << num_edges_1 << std::endl;
        char c;
        unsigned long int st, v;
        bool tab = true;
        if ( tab ) {
            for ( unsigned int i = 0 ; i < num_vertex ; i++ ) {
                ifs >> st >> v;
                new_id[st] = v;
            }
        } else {
            for ( unsigned int i = 0 ; i < num_vertex ; i++ ) {
                ifs >> c >> st >> c >> v >> c;
                new_id[st] = v;
            }
        }
        ifs.close();
    }

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

    void reorder(Algo algo){
        switch(algo) {
            case Algo::hisorder_wo_blc:
                std::cout << "[!!] reordering method: [ history-based order ]" << '\n';
                HisOrder_wo_blc();
                break;
            case Algo::hisorder:
                std::cout << "[!!] reordering method: [ history-based order + balance ]" << '\n';
                HisOrder();
                break;
            case Algo::hisorder_pcpm:
                std::cout << "[!!] Reordering method: [ history-based order + balance(for pcpm) ]" << std::endl;
                HisOrder_PCPM();
                break;
            case Algo::randm:
                std::cout << "reordering method: random" << '\n';
                fastRandom();
                break;
            case Algo::sort:
                std::cout << "reordering method: sort" << '\n';
                fastSort();
                break;
            case Algo::fbc:
                std::cout << "reordering method: fbc" << '\n';
                fastFBC();
                break;
            case Algo::hc:
                std::cout << "reordering method: hc" << '\n';
                fastHC();
                break;
            case Algo::dbg:
                std::cout << "reordering method: dbg" << '\n';
                fastDBG(8);
                break;
            case Algo::corder:
                std::cout << "reordering method: corder" << '\n';
                fastCorder();
                break;
            case Algo::map:
                std::cout << "reorder according to mapping file(rabbit and gorder)" << '\n';
                mapReorder(this->graph->in_feat);
                break;
            default:
                std::cout << "choose a correct algorithm!" << '\n';
        }
    }
};