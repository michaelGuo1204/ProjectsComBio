# ProjectsComBio
Repo for Computational biology Spring 2024

If you want a neighboring graph with curvature computed, use 

```python
from model.neighbors import neighbors
neighbors(ad_data,n_neighbors=30,n_initial_neighbors=30,...,)
```

If you want to trim the neighbor graph, just select a smaller n_neighbors like 

```python
from model.neighbors import neighbors
neighbors(ad_data,n_neighbors=10,n_initial_neighbors=30,...,)
```

All reproducible notebooks are placed in `Experiment` directory