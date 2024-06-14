from unittest import TestCase
from model.ricci_curv import ricciCurvature, groupSinkorn
import numpy as np
import ot
import scipy.sparse as sp
from sklearn.neighbors import KNeighborsTransformer
class Test(TestCase):

    def test_group_sinkorn(self):
        try :
            import cupy as nx
        except ImportError:
            import numpy as nx
        a = nx.array([[.5], [.5]])
        b = nx.array([[.5,.5],[.5,.5]])
        C = nx.array([[0,1],[1,0]],dtype=nx.float32)
        true = ot.sinkhorn2(a,b,C,1e-1)
        res,isvalid = groupSinkorn(a,b,C,1e-1)
        self.assertTrue(nx.allclose(res,true))
    def test_ricci_curvature(self):
        node_feature = np.random.randn(500,100)
        #node_feature = node_feature/np.linalg.norm(node_feature,axis=1).reshape(-1,1)
        graph_dist = KNeighborsTransformer(n_neighbors=15).fit_transform(node_feature)
        graph_index = graph_dist.tocoo().coords[1].reshape(-1,16)
        res = ricciCurvature(graph_dist,graph_index)
        self.assertFalse(np.isnan(res.min()) or np.isnan(res.max()))
