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
import os
os.environ["OMP_NUM_THREADS"] = "1"
import warnings
from sklearn.exceptions import ConvergenceWarning

# sklearn
# from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
# from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, silhouette_samples
from joblib import Parallel, delayed
# local
from core.cluster import Cluster
import core.tools as tools
os.environ["OMP_NUM_THREADS"] = "12"

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



# Отключаем предупреждения ConvergenceWarning
warnings.simplefilter("ignore", category=ConvergenceWarning)

def perform_clustering_n_random(data, krange):
    """Кластеризация с перебором значений K"""
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)

    min_k, max_k, step_k = krange
    set_k = list(range(min_k, max_k + 1, step_k))
    clusters_data = []

    for k in set_k:
        if k > len(data):
            continue
        kmeans = MiniBatchKMeans(n_clusters=k, random_state=42, batch_size=100).fit(data)
        clusters_data.append(kmeans)

    return clusters_data, set_k

def extract_best_k(clusters_data, set_k):
    if not clusters_data:
        raise ValueError("extract_best_k: clusters_data пустой!")

    k_scores = [{"k": set_k[i], "mean": clus.inertia_} for i, clus in enumerate(clusters_data)]
    best_k = min(k_scores, key=lambda x: x["mean"])['k']
    return best_k, k_scores

def pick_best_candidate(clusters_data, best_k):
    for clus in clusters_data:
        if clus.n_clusters == best_k:
            return clus
    raise ValueError(f"Не найдено подходящих кластеров для best_k={best_k}")

def get_representatives(params, coulomb, energies, pfiles):
    """Оптимизированная версия с ускоренным выбором представителей"""
    if np.all(energies == energies[0]):  
        energies = np.array([n if n is not None else sys.maxsize for n in energies])
    if len(params["KMEANS"]) == 1:
        k = min(params["KMEANS"][0], len(coulomb))
        kmeans = MiniBatchKMeans(n_clusters=max(k, 2), random_state=42, batch_size=100).fit(coulomb)
    else:
        clusters_data, set_k = perform_clustering_n_random(np.array(coulomb), params["KMEANS"])
        best_k, _ = extract_best_k(clusters_data, set_k)
        kmeans = pick_best_candidate(clusters_data, best_k)

    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_

    # Параллельный выбор ближайших точек к центрам кластеров
    def find_closest(i):
        cluster_indices = np.where(labels == i)[0]
        if len(cluster_indices) == 0:
            return None
        distances = np.linalg.norm(coulomb[cluster_indices] - cluster_centers[i], axis=1)
        return cluster_indices[np.argmin(distances)]

    representative_indices = Parallel(n_jobs=-1)(delayed(find_closest)(i) for i in range(len(cluster_centers)))
    return [idx for idx in representative_indices if idx is not None]

