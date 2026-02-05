import pandas as pd
from tcr_seq_analysis.profile_plot import ProfilePlot
from test.setup import TestCase


class TestProfilePlot(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()
    
    def test_main(self):
        ProfilePlot(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count_df.csv', index_col=0)
        )
