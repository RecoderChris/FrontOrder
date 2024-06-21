# clustering dataset
# determine k using elbow method

from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import numpy as np
import matplotlib.pyplot as plt
import argparse

def read_feature_vectors(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    feature_vectors = []
    for line in lines:
        vector = [int(x) for x in line.strip().split(',')]
        feature_vectors.append(vector)
    return np.array(feature_vectors)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, default="R20", help="dataset name")
    parser.add_argument("--feat_file", type=str, default="../dataset/feature/R20_dim_10.feat", help="feature dimension(<=the number of sampling)")
    args = parser.parse_args()
    
    file_path = args.feat_file
    data = args.data
    print("file path = {}".format(file_path))
    
    X = read_feature_vectors(file_path)

    print(f"==== {data} ====")

    plt.plot()
    plt.title('Dataset')
    print(X.shape)
    colors = ['b', 'g', 'r']
    markers = ['o', 'v', 's']

    # k means determine k
    distortions = []
    K = np.logspace(2, 6, 5, base=2).astype(int)
    for k in K:
        print(f"calculate {k}")
        kmeanModel = KMeans(n_clusters=k).fit(X)
        distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0])

    # Plot the elbow
    plt.plot(K, distortions, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Distortion')
    plt.title('The Elbow Method showing the optimal k')
    plt.savefig(f"elbow-{data}.png")
    plt.clf()