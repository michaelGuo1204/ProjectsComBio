import scanpy.logging as logg
import numpy as np
from typing import Union
from anndata import AnnData
AnyRandom = Union[int, np.random.RandomState, None]
from .TrimNG import TrimNeighbors,_Trim_Scheme
# interface to scanpy.pp
def neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    n_initial_neighbors: int | None = None,
    scheme: _Trim_Scheme|None = None,
    n_pcs: int | None = None,
    *,
    method = "umap",
    use_rep: str | None = None,
    knn: bool = True,
    random_state: int = 0,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | None:

    start = logg.info("computing neighbors")
    adata = adata.copy() if copy else adata
    if adata.is_view:  # we shouldn't need this here...
        adata._init_as_actual(adata.copy())
    neighbors = TrimNeighbors(adata)
    neighbors.compute_neighbors(
        n_neighbors=n_neighbors,
        n_neighbors_set=n_initial_neighbors,
        trimming_scheme=scheme,
        n_pcs=n_pcs,
        use_rep=use_rep,
        knn=knn,
        random_state=random_state,
    )

    if key_added is None:
        key_added = "neighbors"
        conns_key = "connectivities"
        dists_key = "distances"
        curvature_key = "curvature"
    else:
        conns_key = key_added + "_connectivities"
        dists_key = key_added + "_distances"
        curvature_key = key_added + "_curvature"

    adata.uns[key_added] = {}

    neighbors_dict = adata.uns[key_added]

    neighbors_dict["connectivities_key"] = conns_key
    neighbors_dict["distances_key"] = dists_key
    neighbors_dict["curvature_key"] = curvature_key
    neighbors_dict["params"] = dict(
        n_neighbors=neighbors.n_neighbors,
        method = method,
        trimming_scheme=scheme,
        random_state=random_state,
    )
    if use_rep is not None:
        neighbors_dict["params"]["use_rep"] = use_rep
    if n_pcs is not None:
        neighbors_dict["params"]["n_pcs"] = n_pcs

    adata.obsp[dists_key] = neighbors.distances
    adata.obsp[conns_key] = neighbors.connectivities
    adata.obsp[curvature_key] = neighbors.curvature

    if neighbors.rp_forest is not None:
        neighbors_dict["rp_forest"] = neighbors.rp_forest
    logg.info(
        "    finished",
        time=start,
        deep=(
            f"added to `.uns[{key_added!r}]`\n"
            f"    `.obsp[{dists_key!r}]`, distances for each pair of neighbors\n"
            f"    `.obsp[{conns_key!r}]`, weighted adjacency matrix\n"
            f"    `.obsp[{curvature_key!r}]`, curvature matrix"

        ),
    )
    return adata if copy else None