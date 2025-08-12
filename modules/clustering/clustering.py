import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.cluster import AgglomerativeClustering
from sklearn_extra.cluster import KMedoids
from sklearn.manifold import SpectralEmbedding
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
import umap
import hdbscan
import matplotlib.pyplot as plt
import seaborn as sns
from prolif.utils import to_bitvectors, to_countvectors


def tanimoto_distance_matrix(fps):
    n = len(fps)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            dist = 1 - sim
            dist_matrix[i, j] = dist_matrix[j, i] = dist
    return dist_matrix

def lead_discovery_clustering(bitvecs, n_clusters=5):
    dist_matrix = tanimoto_distance_matrix(bitvecs)
    model = AgglomerativeClustering(n_clusters=n_clusters, affinity="precomputed", linkage="average")
    labels = model.fit_predict(dist_matrix)
    return labels

def diversity_analysis_clustering(bitvecs, n_clusters=10):
    dist_matrix = tanimoto_distance_matrix(bitvecs)
    model = KMedoids(n_clusters=n_clusters, metric="precomputed", init="k-medoids++", random_state=42)
    labels = model.fit_predict(dist_matrix)
    return labels

def hit_expansion_clustering(countvecs):
    reducer = umap.UMAP(metric="euclidean")
    embedding = reducer.fit_transform(countvecs)
    clusterer = hdbscan.HDBSCAN(min_cluster_size=5)
    labels = clusterer.fit_predict(embedding)
    return embedding, labels