from tcr_seq_analysis.compile_table import CompileTable
from test.setup import TestCase


class TestCompileTable(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()
    
    def test_main(self):
        count_df = CompileTable(self.settings).main(
            csv_dir=f'{self.indir}/csv_dir',
            csv_suffix='_TRB.UMI_1.immune_viewer_report.csv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            clonal_index_column='cdr3_amino_acid_sequence',
            count_column='read_count',
        )
        count_df.to_csv(f'{self.outdir}/count_df.csv', index=True)
        self.assertFileEqual(
            f'{self.outdir}/count_df.csv',
            f'{self.indir}/count_df.csv'
        )
