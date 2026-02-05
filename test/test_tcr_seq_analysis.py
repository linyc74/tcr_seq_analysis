from tcr_seq_analysis.tcr_seq_analysis import TcrSeqAnalysis
from test.setup import TestCase


class TestTcrSeqAnalysis(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        TcrSeqAnalysis(self.settings).main(
            csv_dir=f'{self.indir}/csv-dir',
            csv_suffix='_TRB.UMI_1.immune_viewer_report.csv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            clonal_index_column='cdr3_amino_acid_sequence',
            count_column='read_count',
            group_column='Group',
            rpm_cutoff=1000,
            clustering_identity=0.60,
            p_value=0.05,
            colormap='Set1',
            invert_colors=False)
