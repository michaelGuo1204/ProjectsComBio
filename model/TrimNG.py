from types import MappingProxyType
from typing import Literal

import numpy as np
from scanpy.neighbors import Neighbors
from scanpy.neighbors._types import _Method, KnnTransformerLike, _KnownTransformer, _Metric, _MetricFn

from typing import  Any, Mapping
from scipy.sparse import csr_matrix,issparse
from scipy.sparse.csgraph import csgraph_to_masked
from anndata import AnnData
from scanpy._utils import NeighborsView, AnyRandom
from scanpy.logging import debug

from .ricci_curv import ricciCurvature
_Trim_Scheme = Literal["min","max"]
class TrimNeighbors(Neighbors):
    def __init__(self,
                 adata:AnnData,
                 neighbors_key:str|None = None):
        super().__init__(adata)
        self._curvature: np.ndarray|csr_matrix|None = None
        if neighbors_key is None:
            neighbors_key = "neighbors"
        if neighbors_key in adata.uns:
            neighbors = NeighborsView(adata, neighbors_key)
            self.knn = issparse(neighbors["curvature"])
            self._curvature = neighbors["curvature"]

    @property
    def curvature(self) -> np.ndarray | csr_matrix | None:
        return self._curvature

    def compute_neighbors(
        self,
        n_neighbors: int = 15,
        n_neighbors_set: int = 15,
        trimming_scheme: _Trim_Scheme|None = None,
        n_pcs: int | None = None,
        *,
        use_rep: str | None = None,
        knn: bool = True,
        method: _Method | None = "umap",
        transformer: KnnTransformerLike | _KnownTransformer | None = None,
        metric: _Metric | _MetricFn = "euclidean",
        metric_kwds: Mapping[str, Any] = MappingProxyType({}),
        random_state: AnyRandom = 0,
    ) -> None:
        super().compute_neighbors(n_neighbors=n_neighbors_set+1, n_pcs=n_pcs, use_rep=use_rep, knn=knn, method=method,
                                  transformer=transformer, metric=metric, metric_kwds=metric_kwds, random_state=random_state)
        start_curvature = debug("computing curvature")
        if self.knn:
            self._curvature = self._compute_curvature(computed_neighbors=n_neighbors_set)
        start_trimming = debug("trimming neighbors")
        if n_neighbors != n_neighbors_set:
            assert trimming_scheme is not None
            self._trimming(trimming_scheme,n_neighbors)


    def _compute_curvature(self,computed_neighbors:int):
        graph_index = self._distances.tocoo().coords[1].reshape(-1,computed_neighbors)
        curvature = ricciCurvature(self._distances,graph_index)
        return curvature

    def _trimming(self,scheme,n_neighbors):
        n_samples = self._distances.shape[0]
        graph_mask = csgraph_to_masked(self._distances)
        curvature_masked = graph_mask * self._curvature
        curvature_graph = csr_matrix((n_samples,n_samples))
        connectivity_graph = csr_matrix((n_samples,n_samples))
        distance_graph = csr_matrix((n_samples,n_samples))
        if scheme == 'max':
            target_curv_indices = np.argsort(np.abs(curvature_masked), axis=1)[:,:n_neighbors]
        elif scheme == 'min':
            target_curv_indices = np.argsort(np.abs(curvature_masked), axis=1)[:,::-1][:,:n_neighbors]
        else:
            raise ValueError(f"Unknown trimming scheme {scheme}")
        for i in range(n_samples):
            curvature_graph[i, target_curv_indices[i]] = self._curvature[i, target_curv_indices[i]]
            connectivity_graph[i, target_curv_indices[i]] = self._connectivities[i, target_curv_indices[i]]
            distance_graph[i, target_curv_indices[i]] = self._distances[i, target_curv_indices[i]]
        self._curvature = curvature_graph
        self._connectivities = connectivity_graph
        self._distances = distance_graph
        self.n_neighbors = n_neighbors
        pass







