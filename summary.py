import json
import numpy as np


data_list = ["DL", "KR", "MPI", "TW", "FR", "RM"] #
reorders=["OG", "SO", "HC", "DBG", "FBC", "CO", "RBT", "GO", "FO"]
algos=["BFS", "SSSP", "WCC", "PageRank"]

directory="result"
def summary_ligra():
    data_folder = f"{directory}/ligra" 
    print("Figure 8.(b) Avg speedups to OG on Ligra: ")
    ligra_dict = {}
    algo_map = {"BFS": "BFS", 
                "SSSP": "BellmanFord", 
                "WCC": "Components", 
                "PageRank" : "PageRankDelta"}
    for key, algo in algo_map.items():
        print("|===============|\n|{}|\n|---------------|".format(key.center(15)))
        ligra_dict[algo] = {}
        for data in data_list:
            # print("-----\ndataset = {}".format(data))
            ligra_dict[algo][data] = {}
            for rod in reorders:
                log_path = f"{data_folder}/{data}/{algo}/{rod}.log"
                latency = 0
                repeat = 0
                try:
                    with open(log_path, "r") as log_file:
                        line = log_file.readline()
                        while line:
                            if "Running" in line.split():
                                latency += float(line.split()[-1])
                                repeat += 1
                            line = log_file.readline()
                    ligra_dict[algo][data][rod] = latency / repeat
                except FileNotFoundError:
                    ligra_dict[algo][data][rod] = 0
                # print("{} = {:.2f}".format(rod, ligra_dict[algo][data][rod]))
        ligra_dict[algo]["AVG"] = {}
        print("|{}|{}|".format("Reorder", "Speedup"))
        for rod in reorders:
            speedup = 0
            cnt = 0
            for data in data_list:
                if(ligra_dict[algo][data][rod] == 0):
                    speedup += 0
                else:
                    speedup += (ligra_dict[algo][data]["OG"] / ligra_dict[algo][data][rod])
                    cnt += 1
            ligra_dict[algo]["AVG"][rod] = speedup / cnt
            float_value = "{:.2f}".format(ligra_dict[algo]["AVG"][rod])
            print("|{}|{}|".format(rod.center(7), float_value.center(7)))
    json_path = "ligra.json".format(data_folder)
    with open(json_path, "w") as f:
        json.dump(ligra_dict, f, indent=4, sort_keys=False)

def summary_gpop():
    print("Figure 8.(a) Avg speedups to OG on GPOP: ")
    data_folder = f"{directory}/GPOP" 
    ligra_dict = {}
    algo_map = {"BFS": "bfs", 
                "SSSP": "sssp", 
                "WCC": "cc", 
                "PageRank" : "pr"}
    for key, algo in algo_map.items():
        print("|===============|\n|{}|\n|---------------|".format(key.center(15)))
        ligra_dict[algo] = {}
        for data in data_list:
            # print("-----\ndataset = {}".format(data))
            ligra_dict[algo][data] = {}
            for rod in reorders:
                log_path = f"{data_folder}/{data}/{algo}/{rod}.log"
                latency = 0
                repeat = 0
                try:
                    with open(log_path, "r") as log_file:
                        line = log_file.readline()
                        while line:
                            if "Avg" in line.split():
                                latency += float(line.split()[-1])
                                repeat += 1
                            line = log_file.readline()
                    ligra_dict[algo][data][rod] = latency / repeat
                except FileNotFoundError:
                    ligra_dict[algo][data][rod] = 0
                # print("{} = {:.2f}".format(rod, ligra_dict[algo][data][rod]))
        ligra_dict[algo]["AVG"] = {}
        print("|{}|{}|".format("Reorder", "Speedup"))
        for rod in reorders:
            speedup = 0
            cnt = 0
            for data in data_list:
                if(ligra_dict[algo][data][rod] == 0):
                    speedup += 0
                else:
                    speedup += (ligra_dict[algo][data]["OG"] / ligra_dict[algo][data][rod])
                    cnt+=1
            ligra_dict[algo]["AVG"][rod] = speedup / cnt
            float_value = "{:.2f}".format(ligra_dict[algo]["AVG"][rod])
            print("|{}|{}|".format(rod.center(7), float_value.center(7)))
    json_path = "gpop.json".format(data_folder)
    with open(json_path, "w") as f:
        json.dump(ligra_dict, f, indent=2, sort_keys=False)


if __name__ == "__main__":
    summary_gpop()
    summary_ligra()