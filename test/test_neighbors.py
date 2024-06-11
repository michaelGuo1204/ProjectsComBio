from unittest import TestCase
import numpy as np
import anndata as ad
import scanpy as sc
from model.neighbors import neighbors

class Test(TestCase):
    def test_neighbors_nontrim(self):
        X = np.random.randn(100, 30)
        ad_data = ad.AnnData(X)
        try:
            neighbors(ad_data, n_neighbors=15, n_initial_neighbors=15, scheme="min")
        except Exception as e:
            self.fail(e)
    def test_neighbors_trimming(self):
        X = np.random.randn(100, 30)
        ad_data = ad.AnnData(X)
        try:
            neighbors(ad_data, n_neighbors=15, n_initial_neighbors=30, scheme="min")
        except Exception as e:
            self.fail(e)

