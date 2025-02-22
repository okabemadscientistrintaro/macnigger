# from ase import Atoms
# from dscribe.descriptors import CoulombMatrix
import numpy as np
import numpy.typing as npt
import tqdm
# import glob
# import math
# import ast
# import os
import sys


# sklearn
# from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
# from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

# local
from core.cluster import Cluster
import core.tools as tools
import os



import numpy as np
from sklearn.cluster import KMeans
import sys
os.environ["OMP_NUM_THREADS"] = "1"

def calculate_silhouette(data: npt.NDArray[np.float64], labels: npt.NDArray[np.int32], k: int) -> float:
    """
    Compute Silhouette for each sample, take the average for each cluster and
    return (#clusters that surpassed the general Silhouette) / (#clusters)
    """
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)
    
    silh_samples = silhouette_samples(data, labels)
    silh_score = silhouette_score(data, labels)

    csilh = np.array(list(zip(labels, silh_samples)))
    csilh = csilh[csilh[:, 0].argsort()]

    groups = np.split(csilh[:, 1], np.unique(csilh[:, 0], return_index=True)[1][1:])

    return sum([group.mean() >= silh_score for group in groups]) / k

def perform_clustering_n_random(data, krange):
    """
    Выполняет кластеризацию данных с разными значениями K.
    """
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)
    
    min_k, max_k, step_k = krange
    set_k = list(range(min_k, max_k + 1, step_k))
    clusters_data = []
    
    for k in set_k:
        if k > len(data):
            continue
        kmeans = KMeans(n_clusters=k, random_state=0, n_init=10).fit(data)
        clusters_data.append(kmeans)
    
    return clusters_data, set_k

def extract_best_k(clusters_data, set_k):
    if not clusters_data:
        raise ValueError("extract_best_k: clusters_data пустой!")
    
    k_scores = []
    for i, clus in enumerate(clusters_data):
        scores = clus.inertia_
        k_scores.append({"k": set_k[i], "mean": scores})
    
    best_k = min(k_scores, key=lambda x: x["mean"])['k']
    return best_k, k_scores

def pick_best_candidate(clusters_data, best_k):
    for clus in clusters_data:
        if clus.n_clusters == best_k:
            return clus
    raise ValueError(f"Не найдено подходящих кластеров для best_k={best_k}")

def get_representatives(params: dict, coulomb, energies, pfiles: list):
    if len(set(energies)) == 1:
        energies = np.array([n if n is not None else sys.maxsize for n in energies])
    
    if len(params["KMEANS"]) == 1:
        k = params["KMEANS"][0]
        k = min(k, len(coulomb))
        kmeans = KMeans(n_clusters=k, random_state=0, n_init=10).fit(coulomb)
    else:
        clusters_data, set_k = perform_clustering_n_random(coulomb, params["KMEANS"])
        best_k, _ = extract_best_k(clusters_data, set_k)
        kmeans = pick_best_candidate(clusters_data, best_k)
    
    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_

    # Выбираем ближайшие точки к центрам кластеров
    representative_indices = []
    for i in range(len(cluster_centers)):
        cluster_indices = np.where(labels == i)[0]  # Индексы точек в кластере i
        if len(cluster_indices) > 0:
            distances = np.linalg.norm(coulomb[cluster_indices] - cluster_centers[i], axis=1)
            rep_index = cluster_indices[np.argmin(distances)]
            representative_indices.append(rep_index)

    return representative_indices

