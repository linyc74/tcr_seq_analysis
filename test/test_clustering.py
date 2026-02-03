import pandas as pd
from tcr_seq_analysis.clustering import Clustering
from test.setup import TestCase


class TestClustering(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        Clustering(self.settings).main(
            df=pd.read_csv(f'{self.indir}/rpm-table.csv', index_col=0),
            abundance_cutoff=500,
            clustering_identity=0.80,
        )
