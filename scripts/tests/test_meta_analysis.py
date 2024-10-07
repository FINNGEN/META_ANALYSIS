#!/usr/bin/env python3

import unittest
from unittest.mock import MagicMock
import numpy as np
import math
from scipy.stats import chi2
from meta_analysis import het_test, n_meta, inv_var_meta, variance_weight_meta, Study, VariantData, do_meta

class TestMetaAnalysis(unittest.TestCase):

    def assertListAlmostEqual(self, list1, list2, places=5):
        """
        Helper method to assert that two lists of floats are almost equal.
        """
        self.assertEqual(len(list1), len(list2), "Lists have different lengths")
        for i, (a, b) in enumerate(zip(list1, list2)):
            with self.subTest(i=i):
                self.assertAlmostEqual(a, b, places=places)

    def test_het_test_no_heterogeneity(self):
        """
        Test het_test with no heterogeneity (all effect sizes equal to meta effect size).
        Expected p-value should be 1.0.
        """
        effs_sizes = [0.5, 0.5, 0.5]
        weights = [1.0, 1.0, 1.0]
        effs_size_meta = 0.5

        p_value = het_test(effs_sizes, weights, effs_size_meta)
        self.assertAlmostEqual(p_value, 1.0, places=5)

    def test_het_test_with_heterogeneity(self):
        """
        Test het_test with some heterogeneity.
        """
        effs_sizes = [0.5, 0.6, 0.4]
        weights = [1.0, 1.0, 1.0]
        effs_size_meta = 0.5

        effs_sizes_array = np.array(effs_sizes)
        weights_array = np.array(weights)
        eff_dev = weights_array * ((effs_sizes_array - effs_size_meta) ** 2)
        sum_eff_dev = np.sum(eff_dev)
        expected_p = chi2.sf(sum_eff_dev, df=len(effs_sizes)-1)

        p_value = het_test(effs_sizes, weights, effs_size_meta)
        self.assertAlmostEqual(p_value, expected_p, places=5)

    def test_n_meta_single_study(self):
        """
        Test n_meta with a single study.
        """
        # Mock Study and VariantData
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        study.effective_size = 100
        variant.beta = 0.5
        variant.z_score = 1.96

        studies = [(study, variant)]
        result = n_meta(studies)

        expected_beta_meta = 0.5  # sum_betas / sum_weights = 0.5*10 /10 =0.5
        self.assertEqual(result[0], expected_beta_meta)
        self.assertIsNone(result[1])  # SE is None
        self.assertTrue(0 <= result[2] <=1)  # p-value between 0 and 1
        self.assertIsNotNone(result[3])  # log10 p-value
        self.assertEqual(result[4], [0.5])  # effs_size_org
        self.assertEqual(result[5], [10.0])  # weights

    def test_inv_var_meta_no_se(self):
        """
        Test inv_var_meta with a study that has no standard error.
        Should return None.
        """
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        variant.se = None
        variant.beta = 0.5
        study.name = "Study1"

        studies = [(study, variant)]
        result = inv_var_meta(studies)
        self.assertIsNone(result)

    def test_inv_var_meta_single_study(self):
        """
        Test inv_var_meta with a single study having a valid standard error.
        """
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        variant.se = 0.1
        variant.beta = 0.5
        study.name = "Study1"

        studies = [(study, variant)]
        result = inv_var_meta(studies)

        expected_inv_var = 1 / (0.1 ** 2)  # 100
        expected_beta_meta = variant.beta  # 0.5
        expected_se_meta = np.sqrt(1 / expected_inv_var)  # 0.1

        self.assertAlmostEqual(result[0], expected_beta_meta, places=5)
        self.assertAlmostEqual(result[1], expected_se_meta, places=5)
        self.assertTrue(0 <= result[2] <= 1)  # p-value between 0 and 1
        self.assertIsNotNone(result[3])  # log10 p-value
        self.assertEqual(result[4], [0.5])  # effs_size_org

        # Use the helper method for weights
        self.assertListAlmostEqual(result[5], [100.0], places=5)  # weights

    def test_variance_weight_meta_no_se(self):
        """
        Test variance_weight_meta with a study that has no standard error.
        Should return None.
        """
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        variant.se = None
        variant.beta = 0.5

        studies = [(study, variant)]
        result = variance_weight_meta(studies)
        self.assertIsNone(result)

    def test_variance_weight_meta_multiple_studies(self):
        """
        Test variance_weight_meta with multiple studies.
        """
        study1 = MagicMock(spec=Study)
        variant1 = MagicMock(spec=VariantData)
        variant1.se = 0.1
        variant1.beta = 0.5
        variant1.z_score = 1.96

        study2 = MagicMock(spec=Study)
        variant2 = MagicMock(spec=VariantData)
        variant2.se = 0.2
        variant2.beta = 0.3
        variant2.z_score = 1.96

        studies = [(study1, variant1), (study2, variant2)]
        result = variance_weight_meta(studies)

        inv_var1 = (1 / 0.1) * variant1.z_score  # 10 * 1.96 = 19.6
        inv_var2 = (1 / 0.2) * variant2.z_score  # 5 * 1.96 = 9.8
        sum_weights = inv_var1 + inv_var2  # 29.4
        expected_beta_meta = (inv_var1 * variant1.beta + inv_var2 * variant2.beta) / sum_weights
        expected_beta_meta = (19.6 * 0.5 + 9.8 * 0.3) / 29.4

        self.assertAlmostEqual(result[0], expected_beta_meta, places=5)
        self.assertIsNone(result[1])  # SE is None
        self.assertTrue(0 <= result[2] <=1)  # p-value
        self.assertIsNotNone(result[3])  # log10 p-value
        self.assertEqual(result[4], [0.5, 0.3])  # effs_size_org
        self.assertEqual(result[5], [19.6, 9.8])  # weights

    def test_n_meta_mismatched_lengths(self):
        """
        Test n_meta with mismatched lengths of effect sizes and studies.
        Should return None.
        """
        # Create multiple studies but provide mismatched effect sizes
        study1 = MagicMock(spec=Study)
        variant1 = MagicMock(spec=VariantData)
        study1.effective_size = 100
        variant1.beta = 0.5
        variant1.z_score = 1.96

        study2 = MagicMock(spec=Study)
        variant2 = MagicMock(spec=VariantData)
        study2.effective_size = 200
        variant2.beta = 0.3
        variant2.z_score = 1.96

        # Intentionally remove one effect size
        studies = [(study1, variant1)]  # Only one study instead of two
        result = n_meta(studies)
        self.assertIsNotNone(result)  # Since len(effs_size) == len(studies), should not be None

    def test_variance_weight_meta_zero_se(self):
        """
        Test variance_weight_meta with a study having zero standard error.
        Should return None.
        """
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        variant.se = 0.0  # Zero standard error
        variant.beta = 0.5

        studies = [(study, variant)]
        result = variance_weight_meta(studies)
        self.assertIsNone(result)

if __name__ == '__main__':
    unittest.main()
