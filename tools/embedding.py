"""
This python script is used to generate feature vector based on the sampled historical frontiers. 
Author: Xinmiao Zhang
email: zhangxinmiao20s@ict.ac.cn
"""

import os
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, default="R20", help="dataset name")
    parser.add_argument("--feat_dim", type=int, default=10, help="feature dimension(<=the number of sampling)")
    parser.add_argument("--data_dir", type=str, default="../dataset", help="directory of dataset ($DATA_DIR)")
    
    args = parser.parse_args()
    print("[Generating feature vector]")
    print("----------\nCheck Configs:")
    print("  data = {}".format(args.data))
    print("  feat_dim = {}".format(args.feat_dim))
    print("  data_dir = {}".format(args.data_dir))
    root_dir = args.data_dir
    sample_dir = "{}/sample".format(root_dir)
    feat_dir = "{}/feature".format(root_dir)
    feat_dim = args.feat_dim
    data = args.data
    if not os.path.exists(feat_dir):
        os.makedirs(feat_dir)
    all_data = []
    print("----------\nReading Sampled Frontiers...")
    for k in list(range(feat_dim)):
        round_data = []
        sample_file = "{}/{}.bfs.{}".format(sample_dir, data, str(k))
        with open(sample_file, "r") as file:
            num_vertex = int(file.readline().strip())  
            real_line = 0
            for line in file:
                line = line.strip()  
                if line != "----------":
                    content = line.split()  
                    content = [int(num) for num in content]
                    round_data.append(content)
                    real_line += 1
            print(f"read({k}) BFS history: niters = {real_line}")
        all_data.append(round_data)
    feat = np.zeros((num_vertex, len(all_data)))
    print("----------\nGenerating Embedding...")
    for round, round_data in enumerate(all_data):
        for iter, iter_data in enumerate(round_data):
            for vid in iter_data:
                if feat[vid][round] != 0:
                    print("wrong!")
                    exit(0)
                feat[vid][round] = iter
    feat_file =  "{}/{}_dim_{}.feat".format(feat_dir, data, len(all_data))
    print("----------\nSaving Feature Vector at {}".format(feat_file))
    np.savetxt(feat_file, feat, delimiter=',', fmt='%d')
    print("Success!\n")
    