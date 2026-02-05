import pandas as pd
from tcr_seq_analysis.clustering import Clustering
from test.setup import TestCase


class TestClustering(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_main(self):
        count_df, rpm_df, motif_count_df, motif_rpm_df = Clustering(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count_df.csv', index_col=0),
            rpm_df=pd.read_csv(f'{self.indir}/rpm_df.csv', index_col=0),
            rpm_cutoff=1000,
            clustering_identity=0.80,
        )
