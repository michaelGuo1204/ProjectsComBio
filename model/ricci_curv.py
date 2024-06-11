import cupy as nx
import numpy as np
from einops import rearrange
from sklearn.metrics import pairwise_distances
from scipy.sparse.csgraph import  shortest_path
from tqdm import tqdm
def groupSinkorn(a, b, C, reg=1e-3,log=True):
    '''
    Compute unbalanced Sinkorn distance between measure a and its neighbors b

    :param a: source measure with shape (n,1)
    :param b: dst measure with shape (n,k), k is the neighbors of a
    :param C: pairwise cost matrix with shape (n,n)
    :param reg: regularity parameter
    :param log: whether to print log file
    :return:
    '''
    assert a.shape[0] == b.shape[0] and a.shape[0] == C.shape[1]
    K = nx.exp(C / (-reg))
    #K -= np.identity(K.shape[0])
    Kp = (1/a).reshape(-1,1) * K
    u, v = nx.ones_like(a) / a.shape[0], nx.ones_like(b) / b.shape[0]
    for ii in range(1,100):
        u_prev = u
        v_prev = v
        KtU = K.T @ u
        v = b / KtU
        u = 1. / (Kp @ v)
        if (nx.any(KtU == 0)
                or nx.any(nx.isnan(u)) or nx.any(nx.isnan(v))
                or nx.any(nx.isinf(u)) or nx.any(nx.isinf(v))):
            # Reached the machine precision
            print('Warning: numerical errors at iteration %d' % ii)
            u = u_prev
            v = v_prev
            break
        if ii % 10 == 0:
            err = (nx.sum((u - u_prev) ** 2) + nx.sum((v - v_prev) ** 2))/u.shape[0]
            if err < 1e-6:
                if log:print("converge")
                break
    res = nx.einsum('ik,ij,jk,ij->k', u, K, v,C)
    return res.get()

def ricciCurvature(knn_dist,knn_index):
    '''
    Compute discrete Ricci Curvature

    :param node_features: np.ndarray for node features
    :param adjacency: Sparse adjacency matrix computed by sklearn or scanpy
    :return: Curvature matrix with same shape as adjacency
    '''
    knn_dist = knn_dist.toarray()
    rw_measure = 1 / knn_dist
    rw_measure = np.nan_to_num(rw_measure, posinf=1e-7)
    rw_measure = (rw_measure / rw_measure.sum(axis=1,keepdims=True))
    rw_measure = nx.asarray(rw_measure)
    cost_matrix = nx.asarray(shortest_path(knn_dist,directed=False))
    curvature = knn_dist.copy()
    samples = knn_dist.shape[0]
    # Parallel computing
    for i in tqdm(range(0, samples)):
        neighbors = knn_index[i]
        measure_i = rw_measure[i].reshape(-1,1)
        dist_neigh = knn_dist[i][neighbors]
        measure_neigh = rw_measure[neighbors].T
        wasserstein_dist = groupSinkorn(a=measure_i, b=measure_neigh, C=cost_matrix, reg=1,log=False)
        curvature_chunk = 1 - wasserstein_dist / (1e-7 + dist_neigh)
        curvature[i][neighbors] = curvature_chunk
    return curvature
