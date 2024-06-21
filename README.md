# Hisorder
具体可以参考`run_hisorder.sh`
```shell
./hisorder -d /path/to/your/data -a 2 -s 1024 -o /path/to/store/data
```

    ./hisorder
    -d [ --data ] arg           input data path
    -o [ --output ] arg         output file path
    -s [ --size ] arg (=1024)   partition size
    -a [ --algorithm ] arg (=0) reordering algorithm
    -t [ --thread ] arg (=20)   threads
    -v [ --vertex ] arg (=100)  start vertex id(for new start vertex id)
    -f [ --feat ] arg (=10)     feature size, for hisorder algorithm
    -k [ --kv ] arg             k value for kmeans
    -i [ --input_feat ] arg     input feature file
    -r [ --output_map]  arg     output mapping file
    注意修改: 目前, -s 已经直接被修改为节点的个数, 而不是节点的KB数. 例如现在-s 1024代表一个分块中含有1024个节点(10.7)
    [ --algorithm ]
    original = 0,
    hisorder_wo_blc = 1,
    hisorder = 2,

