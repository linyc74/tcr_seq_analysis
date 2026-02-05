import pandas as pd
from tcr_seq_analysis.differential_abundance import DifferentialAbundance
from test.setup import TestCase


class TestDifferentialAbundance(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        DifferentialAbundance(self.settings).main(
            df=pd.read_csv(f'{self.indir}/motif_count_df.csv', index_col=0),
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            group_column='Group',
            colors=[(0, 0, 1, 1), (1, 0, 0, 1), (0, 1, 0, 1)],
            p_value=0.05,
        )
