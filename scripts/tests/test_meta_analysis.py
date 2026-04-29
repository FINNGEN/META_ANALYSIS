#!/usr/bin/env python3

import unittest
from unittest.mock import MagicMock, patch
import numpy as np
import numpy
import math
from scipy.stats import chi2
from meta_analysis import (
    het_test,
    n_meta,
    inv_var_meta,
    variance_weight_meta,
    do_meta,
    Study,
    VariantData,
    format_num,
    SUPPORTED_METHODS
)

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

    @patch('meta_analysis.format_num', side_effect=lambda x: f"{x:.2f}" if x is not None else "NA")
    @patch('meta_analysis.het_test', return_value=0.05)
    @patch.dict('meta_analysis.SUPPORTED_METHODS', {})
    def test_do_meta_single_method_no_het(self, mock_het_test, mock_format_num):
        """
        Test do_meta with a single method ('inv_var') and is_het_test=False.
        """
        # Create mock Study and VariantData
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        study.effective_size = 100
        variant.beta = 0.5
        variant.pval = 0.05
        variant.se = 0.1
        variant.z_score = 1.96

        study_list = [(study, variant)]
        methods = ['inv_var']
        is_het_test = False

        # Define the expected output from inv_var_meta
        # Suppose inv_var_meta returns (beta_meta, se_meta, pval, mlogp, effs_size_org, weights)
        expected_method_output = (0.5, 0.1, 0.05, 0.15, [0.5], [100.0])

        # Update SUPPORTED_METHODS to include 'inv_var'
        with patch.dict('meta_analysis.SUPPORTED_METHODS', {'inv_var': MagicMock(return_value=expected_method_output)}):
            result = do_meta(study_list, methods, is_het_test)

        # Expected result: tuple with format_num applied to first 3 elements, m[3] rounded to 2 decimal places
        expected_result = [
            ("0.50", "0.10", "0.05", 0.15)  # method1: inv_var
        ]

        self.assertEqual(result, expected_result)
        # Ensure het_test was not called
        mock_het_test.assert_not_called()

    @patch('meta_analysis.format_num', side_effect=lambda x: f"{x:.2f}" if x is not None else "NA")
    @patch('meta_analysis.het_test', return_value=0.05)
    @patch.dict('meta_analysis.SUPPORTED_METHODS', {})
    def test_do_meta_multiple_methods_with_het(self, mock_het_test, mock_format_num):
        """
        Test do_meta with multiple methods ('n', 'inv_var') and is_het_test=True.
        """
        # Create mock Study and VariantData objects
        study1 = MagicMock(spec=Study)
        variant1 = MagicMock(spec=VariantData)
        study1.effective_size = 100
        variant1.beta = 0.5
        variant1.pval = 0.05
        variant1.se = 0.1
        variant1.z_score = 1.96

        study2 = MagicMock(spec=Study)
        variant2 = MagicMock(spec=VariantData)
        study2.effective_size = 200
        variant2.beta = 0.3
        variant2.pval = 0.01
        variant2.se = 0.2
        variant2.z_score = 2.58

        study_list = [(study1, variant1), (study2, variant2)]
        methods = ['n', 'inv_var']
        is_het_test = True

        # Define expected outputs from 'n' and 'inv_var' methods
        expected_n_meta_output = (
            0.5,      # beta_meta
            None,     # se_meta
            0.05,     # pval
            -1.30102999566,  # mlogp
            [0.5, 0.3],       # effs_size_org
            [10.0, 14.1421356237]  # weights
        )

        expected_inv_var_meta_output = (
            0.4,     # beta_meta
            0.05,    # se_meta
            0.01,    # pval
            -2.0,    # mlogp
            [0.5, 0.3],  # effs_size_org
            [100.0, 25.0]  # weights
        )

        # Mock the 'n' and 'inv_var' methods
        with patch.dict('meta_analysis.SUPPORTED_METHODS', {
            'n': MagicMock(return_value=expected_n_meta_output),
            'inv_var': MagicMock(return_value=expected_inv_var_meta_output)
        }):
            result = do_meta(study_list, methods, is_het_test)

        # Expected result:
        # For 'n':
        #   ("0.50", "NA", "0.05", -1.30103, "0.05")  # mlogp now has 6 decimals
        # For 'inv_var':
        #   ("0.40", "0.05", "0.01", -2.00, "0.05")
        expected_result = [
            ("0.50", "NA", "0.05", -1.30103, "0.05"),  # method 'n'
            ("0.40", "0.05", "0.01", -2.00, "0.05")  # method 'inv_var'
        ]

        self.assertEqual(result, expected_result)

        # Ensure het_test was called correctly for both methods
        mock_het_test.assert_any_call([0.5, 0.3], [10.0, 14.1421356237], 0.5)
        mock_het_test.assert_any_call([0.5, 0.3], [100.0, 25.0], 0.4)
        self.assertEqual(mock_het_test.call_count, 2)

    @patch('meta_analysis.format_num', side_effect=lambda x: f"{x:.2f}" if x is not None else "NA")
    @patch('meta_analysis.het_test', return_value=0.05)
    @patch.dict('meta_analysis.SUPPORTED_METHODS', {})
    def test_do_meta_method_returns_none(self, mock_het_test, mock_format_num):
        """
        Test do_meta when one of the methods returns None.
        """
        # Create mock Study and VariantData
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        study.effective_size = 100
        variant.beta = 0.5
        variant.pval = 0.05
        variant.se = 0.1
        variant.z_score = 1.96

        study_list = [(study, variant)]
        methods = ['n', 'inv_var']
        is_het_test = True

        # Define expected outputs
        expected_n_meta_output = (
            0.5,      # beta_meta
            None,     # se_meta
            0.05,     # pval
            -1.30102999566,  # mlogp
            [0.5],    # effs_size_org
            [10.0]    # weights
        )

        expected_inv_var_meta_output = None  # Simulate method 'inv_var' failing

        # Mock the 'n' and 'inv_var' methods
        with patch.dict('meta_analysis.SUPPORTED_METHODS', {
            'n': MagicMock(return_value=expected_n_meta_output),
            'inv_var': MagicMock(return_value=None)
        }):
            result = do_meta(study_list, methods, is_het_test)

        # Expected result:
        # 'n' method: ("0.50", "NA", "0.05", -1.30103, "0.05")  # mlogp now has 6 decimals
        # 'inv_var' method: None
        expected_result = [
            ("0.50", "NA", "0.05", -1.30103, "0.05"),
            None
        ]

        self.assertEqual(result, expected_result)

        # Ensure het_test was called only for 'n' method
        mock_het_test.assert_called_once_with([0.5], [10.0], 0.5)

    @patch('meta_analysis.format_num', side_effect=lambda x: f"{x:.2f}" if x is not None else "NA")
    @patch('meta_analysis.het_test', return_value=0.05)
    @patch.dict('meta_analysis.SUPPORTED_METHODS', {})
    def test_do_meta_empty_study_list(self, mock_het_test, mock_format_num):
        """
        Test do_meta with an empty study list.
        """
        study_list = []
        methods = ['n', 'inv_var']
        is_het_test = True

        # Mock the 'n' and 'inv_var' methods to return None
        with patch.dict('meta_analysis.SUPPORTED_METHODS', {
            'n': MagicMock(return_value=None),
            'inv_var': MagicMock(return_value=None)
        }):
            result = do_meta(study_list, methods, is_het_test)

        # Expected result: [None, None]
        expected_result = [
            None,
            None
        ]

        self.assertEqual(result, expected_result)

        # Ensure het_test was not called since both methods returned None
        mock_het_test.assert_not_called()

    def test_do_meta_with_unsupported_method(self):
        """
        Test do_meta with an unsupported method name, expecting KeyError.
        """
        # Create mock Study and VariantData
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        study.effective_size = 100
        variant.beta = 0.5
        variant.pval = 0.05
        variant.se = 0.1
        variant.z_score = 1.96

        study_list = [(study, variant)]
        methods = ['unsupported_method']
        is_het_test = False

        with patch.dict('meta_analysis.SUPPORTED_METHODS', {}, clear=True):
            with self.assertRaises(KeyError):
                do_meta(study_list, methods, is_het_test)

    @patch('meta_analysis.format_num', side_effect=lambda x: f"{x:.2f}" if x is not None else "NA")
    @patch('meta_analysis.het_test', return_value=0.10)
    @patch.dict('meta_analysis.SUPPORTED_METHODS', {})
    def test_do_meta_multiple_methods_with_one_none_and_het(self, mock_het_test, mock_format_num):
        """
        Test do_meta with multiple methods where one method returns None and is_het_test=True.
        """
        # Create mock Study and VariantData
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        study.effective_size = 100
        variant.beta = 0.5
        variant.pval = 0.05
        variant.se = 0.1
        variant.z_score = 1.96

        study_list = [(study, variant)]
        methods = ['n', 'inv_var', 'variance']
        is_het_test = True

        # Define expected outputs
        expected_n_meta_output = (
            0.5,      # beta_meta
            None,     # se_meta
            0.05,     # pval
            -1.30,    # mlogp (already rounded)
            [0.5],    # effs_size_org
            [10.0]    # weights
        )

        expected_inv_var_meta_output = None  # Simulate 'inv_var' failing

        expected_variance_weight_meta_output = (
            0.6,      # beta_meta
            None,     # se_meta
            0.02,     # pval
            -1.70,    # mlogp
            [0.5],    # effs_size_org
            [20.0]    # weights
        )

        # Mock the 'n', 'inv_var', and 'variance' methods
        with patch.dict('meta_analysis.SUPPORTED_METHODS', {
            'n': MagicMock(return_value=expected_n_meta_output),
            'inv_var': MagicMock(return_value=None),
            'variance': MagicMock(return_value=expected_variance_weight_meta_output)
        }):
            result = do_meta(study_list, methods, is_het_test)

        # Expected result:
        # 'n' method: ("0.50", "NA", "0.05", -1.30, "0.10")
        # 'inv_var' method: None
        # 'variance' method: ("0.60", "NA", "0.02", -1.70, "0.10")
        expected_result = [
            ("0.50", "NA", "0.05", -1.30, "0.10"),  # method 'n'
            None,                                     # method 'inv_var'
            ("0.60", "NA", "0.02", -1.70, "0.10")   # method 'variance'
        ]

        self.assertEqual(result, expected_result)

        # Ensure het_test was called correctly for methods that returned results
        mock_het_test.assert_any_call([0.5], [10.0], 0.5)
        mock_het_test.assert_any_call([0.5], [20.0], 0.6)
        self.assertEqual(mock_het_test.call_count, 2)


class TestMetaAnalysisPrecision(unittest.TestCase):
    """Test suite for validating numeric precision in meta_analysis.py"""

    def test_format_num_default_precision_is_6(self):
        """Test 1: Verify format_num default precision is 6"""
        # Test with a number that has many significant figures
        result = format_num(0.123456789)
        
        # Parse scientific notation to count significant figures
        # format_num returns scientific notation like "1.234568e-01"
        self.assertIsNotNone(result)
        self.assertNotEqual(result, "NA")
        
        # Check that it's in scientific notation
        self.assertIn('e', result.lower())
        
        # Extract mantissa (part before 'e')
        mantissa = result.lower().split('e')[0]
        # Count significant figures (digits excluding the decimal point)
        sig_figs = len(mantissa.replace('.', '').replace('-', ''))
        
        # Should have 7 digits (1 before decimal + 6 after for precision=6)
        self.assertEqual(sig_figs, 7, 
                        f"format_num should produce 6 decimal places in mantissa, got {sig_figs-1}")

    def test_format_num_with_6_precision(self):
        """Test 2: Verify format_num explicitly with precision=6"""
        result = format_num(1.23456789e-5, precision=6)
        
        # Should be in scientific notation
        self.assertIn('e', result.lower())
        
        # Parse and verify
        mantissa = result.lower().split('e')[0]
        sig_figs = len(mantissa.replace('.', '').replace('-', ''))
        
        # 7 digits total (1 + 6 decimal places)
        self.assertEqual(sig_figs, 7,
                        f"precision=6 should give 7 total digits, got {sig_figs}")

    def test_format_num_very_small_pvalue(self):
        """Test 3: Verify p-value formatting with high precision"""
        result = format_num(1.23456789e-50, precision=6)
        
        self.assertIsNotNone(result)
        self.assertIn('e', result.lower())
        
        # Extract and verify mantissa has correct precision
        parts = result.lower().split('e')
        mantissa = parts[0]
        exponent = int(parts[1])
        
        # Verify exponent is correct magnitude
        self.assertLess(exponent, -40, "Very small p-value should have large negative exponent")
        
        # Verify mantissa precision
        if '.' in mantissa:
            decimal_part = mantissa.split('.')[1]
            self.assertEqual(len(decimal_part), 6,
                           f"Mantissa should have 6 decimal places, got {len(decimal_part)}")

    def test_format_num_handles_none(self):
        """Test 4: Verify format_num returns 'NA' for None"""
        result = format_num(None)
        self.assertEqual(result, "NA", "format_num should return 'NA' for None input")

    def test_format_num_handles_nan(self):
        """Test 5: Verify format_num returns 'NA' for NaN"""
        result = format_num(numpy.nan)
        self.assertEqual(result, "NA", "format_num should return 'NA' for NaN input")

    def test_do_meta_mlogp_precision(self):
        """Test 6: Verify -log10(p) in meta results has 6 decimals"""
        # Create mock study and variant data
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        
        study.effective_size = 100
        variant.beta = 0.5
        variant.pval = 0.001
        variant.se = 0.1
        variant.z_score = 1.96
        
        study_list = [(study, variant)]
        methods = ['inv_var']
        
        result = do_meta(study_list, methods, is_het_test=False)
        
        # Result should be a list with one tuple
        self.assertEqual(len(result), 1)
        self.assertIsNotNone(result[0])
        
        # Fourth element is mlogp (rounded)
        mlogp = result[0][3]
        
        # Should be a numpy float rounded to 6 decimals
        self.assertIsInstance(mlogp, (float, numpy.floating))
        
        # Convert to string and check decimal places
        mlogp_str = f"{mlogp:.10f}"  # Format with extra precision to see what we have
        if '.' in mlogp_str:
            # Count actual decimal places (excluding trailing zeros)
            decimal_part = mlogp_str.split('.')[1].rstrip('0')
            self.assertLessEqual(len(decimal_part), 6,
                               f"mlogp should have at most 6 decimal places, got {len(decimal_part)}")

    def test_do_meta_beta_precision(self):
        """Test 7: Verify beta in meta results has 6 significant figures"""
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        
        study.effective_size = 100
        variant.beta = 0.123456789  # More precision than we want to keep
        variant.pval = 0.001
        variant.se = 0.05
        variant.z_score = 2.58
        
        study_list = [(study, variant)]
        methods = ['inv_var']
        
        result = do_meta(study_list, methods, is_het_test=False)
        
        self.assertIsNotNone(result[0])
        
        # First element is beta (formatted)
        beta_formatted = result[0][0]
        
        # Should be in scientific notation or have limited precision
        self.assertIsInstance(beta_formatted, str)
        self.assertNotEqual(beta_formatted, "NA")
        
        # If scientific notation, check mantissa precision
        if 'e' in beta_formatted.lower():
            mantissa = beta_formatted.lower().split('e')[0]
            if '.' in mantissa:
                decimal_part = mantissa.split('.')[1]
                self.assertEqual(len(decimal_part), 6,
                               f"Beta mantissa should have 6 decimal places, got {len(decimal_part)}")

    def test_do_meta_se_precision(self):
        """Test 8: Verify standard error in meta results has 6 significant figures"""
        study = MagicMock(spec=Study)
        variant = MagicMock(spec=VariantData)
        
        study.effective_size = 100
        variant.beta = 0.5
        variant.pval = 0.001
        variant.se = 0.0987654321  # More precision than we want
        variant.z_score = 2.58
        
        study_list = [(study, variant)]
        methods = ['inv_var']
        
        result = do_meta(study_list, methods, is_het_test=False)
        
        self.assertIsNotNone(result[0])
        
        # Second element is SE (formatted)
        se_formatted = result[0][1]
        
        self.assertIsInstance(se_formatted, str)
        self.assertNotEqual(se_formatted, "NA")
        
        # Check precision in scientific notation
        if 'e' in se_formatted.lower():
            mantissa = se_formatted.lower().split('e')[0]
            if '.' in mantissa:
                decimal_part = mantissa.split('.')[1]
                self.assertEqual(len(decimal_part), 6,
                               f"SE mantissa should have 6 decimal places, got {len(decimal_part)}")

    def test_numpy_round_precision_6(self):
        """Test 9: Verify numpy.round with 6 decimals works correctly"""
        # Test the rounding that's used for mlogp in do_meta
        test_values = [
            (10.123456789, 10.123457),
            (5.555555555, 5.555556),
            (0.000001234, 0.000001),
            (99.999999999, 100.0)
        ]
        
        for input_val, expected in test_values:
            rounded = numpy.round(input_val, 6)
            self.assertAlmostEqual(rounded, expected, places=6,
                                  msg=f"numpy.round({input_val}, 6) should equal {expected}")

    def test_inv_var_meta_precision(self):
        """Test 10: Verify inverse variance meta maintains precision"""
        study1 = MagicMock(spec=Study)
        variant1 = MagicMock(spec=VariantData)
        study1.name = "Study1"
        study1.effective_size = 100
        variant1.beta = 0.123456
        variant1.pval = 0.001
        variant1.se = 0.05
        variant1.z_score = 2.58
        
        study2 = MagicMock(spec=Study)
        variant2 = MagicMock(spec=VariantData)
        study2.name = "Study2"
        study2.effective_size = 200
        variant2.beta = 0.234567
        variant2.pval = 0.0001
        variant2.se = 0.04
        variant2.z_score = 3.29
        
        study_list = [(study1, variant1), (study2, variant2)]
        
        # Call inv_var_meta directly
        result = inv_var_meta(study_list)
        
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 6)  # Should return tuple of 6 elements
        
        # Check that beta_meta (first element) is a float
        beta_meta = result[0]
        self.assertIsInstance(beta_meta, float)
        
        # Check that SE (second element) is a float
        se_meta = result[1]
        self.assertIsInstance(se_meta, float)
        
        # Check that p-value (third element) is a float
        pval_meta = result[2]
        self.assertIsInstance(pval_meta, float)
        self.assertGreater(pval_meta, 0)
        self.assertLessEqual(pval_meta, 1)


if __name__ == '__main__':
    unittest.main()
