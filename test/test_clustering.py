import pandas as pd
from tcr_seq_analysis.clustering import Clustering
from test.setup import TestCase


class TestClustering(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        count_df, motif_count_df = Clustering(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count_df_in.csv', index_col=0),
            rpm_cutoff=1000,
            clustering_identity=0.80,
        )
        count_df.to_csv(f'{self.outdir}/count_df_out.csv', index=True)
        motif_count_df.to_csv(f'{self.outdir}/motif_count_df.csv', index=True)
        self.assertFileEqual(
            f'{self.outdir}/count_df_out.csv',
            f'{self.indir}/count_df_out.csv'
        )
        self.assertFileEqual(
            f'{self.outdir}/motif_count_df.csv',
            f'{self.indir}/motif_count_df.csv'
        )