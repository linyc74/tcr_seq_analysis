import numpy as np
import pandas as pd
from tcr_seq_analysis.diversity_clonality import DiversityClonality, gini_coefficient, \
    normalized_shannon, observed_clones, shannon_entropy, simpson_concentration
from test.setup import TestCase


class TestDiversityClonality(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        self.settings.for_publication = True
        DiversityClonality(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count_df.csv', index_col=0),
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            group_column='Group',
            colors=[(0, 0, 1, 1), (1, 0, 0, 1), (0, 1, 0, 1)],
        )


class TestFunctions(TestCase):

    def test_gini_coefficient_removes_zeros(self):
        values = np.array([0.0, 2.0, 1.0])
        result = gini_coefficient(values=values)
        self.assertAlmostEqual(result, -1.0 / 3.0)

    def test_observed_clones(self):
        values = np.array([0.0, 1.0, 2.0, 0.0])
        self.assertEqual(observed_clones(values=values), 2)

    def test_shannon_entropy_zero_total(self):
        values = np.array([0.0, 0.0, 0.0])
        self.assertEqual(shannon_entropy(values=values), 0.0)

    def test_shannon_entropy_known_value(self):
        values = np.array([1.0, 1.0])
        expected = float(np.log(2.0))
        self.assertAlmostEqual(shannon_entropy(values=values), expected)

    def test_normalized_shannon_two_clones(self):
        values = np.array([1.0, 1.0])
        self.assertAlmostEqual(normalized_shannon(values=values), 1.0)

    def test_normalized_shannon_single_clone(self):
        values = np.array([5.0, 0.0])
        self.assertEqual(normalized_shannon(values=values), 0.0)

    def test_simpson_concentration(self):
        values = np.array([1.0, 2.0])
        self.assertAlmostEqual(simpson_concentration(values=values), 5/9)
