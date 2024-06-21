# Hisorder
Frontier-guided Distribution Graph Reordering

## What is it?
HisOrder is a graph reordering method, which improves graph locality to reduce the cache misses in graph processing. 

Unlike the previous reordering which depends on static characteristics in graph, HisOrder firstly profiles the traces of **His**torical graph processing (primarily the concurrent frontiers) to construct the locality metric between vertices, and then utilizes an unsupervised ML method (K-means at present) to excavate the clusters of high locality to guide graph **Reorder**ing. Furthermore, since the learned clusters of vertices are more likely to be co-activated, HisOrder also fine-tunes the load balance in parallel graph processing with the clusters. 

![hisorder](img/hisorder.png)

For more details, please refer to [our paper](https://liu-cheng.github.io/). 

## Getting Started
### 0. Dependencies
At the minimum, HisOrder depends on the following software:
- make>=4.1
- g++>=7.5.0 (compliant with C++11 standard and OpenMP is required)
- Boost (>=1.58.0)

A possible command to satisfy these dependencies:
```shell
apt-get install build-essentials libboost-dev libomp-dev
```

### 1. Compilation
To compile HisOrder, run the following on root directory:
```shell
make
```

### 2. Data Preparation
HisOrder takes graph in Compressed Sparse Row (CSR) format as both the input and output. 
In this way, after downloading or generating the graph in edgelist format (i.e. `src_vertex dst_vertex` format), we should convert the data into CSR format firstly for prepration. 
#### 2.1 Download Graph or Generate Graph
First of all, create a directory to preserve all the graph data. 
```shell
ROOT_DIR=`pwd`
TOOL_DIR=$ROOT_DIR/tools
DATA_DIR=$ROOT_DIR/dataset
mkdir -p dataset && cd dataset
```
To obtain the graph data, you can download real-world graph data in edgelist format. 
For instance, to download `soc-LiveJournal` graph from [SNAP large network collection](http://snap.stanford.edu/data/index.html): 
```shell
wget http://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gzip -d soc-LiveJournal1.txt.gz
```
Or you could generate RMAT graph using [`PaRMAT`](https://github.com/farkhor/PaRMAT). To generate a RMAT graph (named `R20.el`) with 1M vertices and 16M edges using 8 threads, you can run: 
```shell
cd ${TOOL_DIR}/PaRMAT/Release
make
./PaRMAT -nEdges 16777216 -nVertices 1048576 -output R20.el -threads 8 -noDuplicateEdges
mv R20.el $DATA_DIR
```
#### 2.2 Convert Edgelist into CSR format
We provide tools to support conversion from edgelist into CSR. 
For example, to convert RMAT graph `R20.el`, you could run:
```shell
cd ${TOOL_DIR}/el2csr
make
./el2csr ${DATA_DIR}/R20.el ${DATA_DIR}/R20.csr
```

### 3. Graph Reordering
```shell
./hisorder -d /path/to/your/data -a 2 -s 1024 -o /path/to/store/data
```

    ./hisorder
    -d [ --data ] arg           input data path
    -o [ --output ] arg         output file path
    -s [ --size ] arg (=1024)   partition size(vertex number)
    -a [ --algorithm ] arg (=0) reordering algorithm
    -t [ --thread ] arg (=20)   threads
    -v [ --vertex ] arg (=100)  start vertex id(for new start vertex id)
    -f [ --feat ] arg (=10)     feature size, for hisorder algorithm
    -k [ --kv ] arg             k value for kmeans
    -i [ --input_feat ] arg     input feature file
    -r [ --output_map]  arg     output mapping file

    [ --algorithm ]
    original = 0,
    hisorder_wo_blc = 1,
    hisorder = 2,

## Evaluation
### Compilation
### Running Graph Algorithms

## Benchmark Summary

## Future Work

## Citation

## References