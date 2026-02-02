from tcr_seq_analysis.compile_table import CompileTable
from test.setup import TestCase


class TestCompileTable(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()
    
    def test_main(self):
        CompileTable(self.settings).main(
            csv_dir=f'{self.indir}/csv_dir',
            csv_suffix='_TRB.UMI_1.immune_viewer_report.csv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            clonal_index_columns=['cdr3_amino_acid_sequence'],
            count_column='read_count',
        )
        self.assertFileEqual(
            f'{self.outdir}/rpm-table.csv',
            f'{self.indir}/rpm-table.csv'
        )
        self.assertFileEqual(
            f'{self.outdir}/count-table.csv',
            f'{self.indir}/count-table.csv'
        )
