# README #

Use the following commands to run

1. make 
2. ./bfs <filename> -s <root node> -t <numThreads(optional)> -rounds <#rounds(default 3)> -map <Addr_of_mappingfile> -data <nameofdata>

./bfs /home/zhangxm/project/data/LJ/LJQ.csr -t 48 -rounds 32 -s 512 -map 

Rounds - repeated runs starting from the root (multiple rounds for timing)
Map - if using *orig* order, one can not specify the map parameter; if using *other* order, one should specify the given address of mapping file. 
 
