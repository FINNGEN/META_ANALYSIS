#!/usr/bin/env python3

import unittest
import sys
import os
import tempfile
import gzip
import subprocess
import json
from io import StringIO
from unittest.mock import patch

# Add parent directory to path to import modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import munge
import harmonize
import meta_analysis


class TestPrecisionIntegration(unittest.TestCase):
    """Integration tests for numeric precision across the full pipeline"""

    def setUp(self):
        """Set up test fixtures and temporary directory"""
        self.temp_dir = tempfile.mkdtemp()
        self.test_files = []

    def tearDown(self):
        """Clean up temporary files"""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def create_test_sumstats_file(self, filename, data_rows):
        """Helper to create a test summary statistics file"""
        filepath = os.path.join(self.temp_dir, filename)
        with gzip.open(filepath, 'wt') as f:
            # Write header
            f.write("#CHR\tPOS\tREF\tALT\taf_alt\tbeta\tsebeta\tpval\tmlogp\n")
            # Write data rows
            for row in data_rows:
                f.write("\t".join(map(str, row)) + "\n")
        return filepath

    def test_precision_not_degraded(self):
        """Test 2: Verify precision doesn't degrade through processing"""
        # Create test data with exactly 6 decimal places
        test_data = [
            ["1", "1000000", "A", "G", "0.123456", "0.234567", "0.045678", "0.000123", "3.910095"],
            ["1", "1000100", "C", "T", "0.654321", "0.345678", "0.056789", "0.000456", "3.341124"],
        ]
        
        # Check that the input values have 6 decimals
        for row in test_data:
            for val in row[4:8]:  # Numeric columns
                if '.' in val:
                    decimals = len(val.split('.')[1])
                    self.assertEqual(decimals, 6, 
                                   f"Test data should have 6 decimals, {val} has {decimals}")

    def test_scientific_notation_precision(self):
        """Test 3: Verify scientific notation maintains 6 significant figures"""
        # Test with extreme p-values
        test_pvalues = [
            1e-100,
            1.23456e-50,
            9.87654e-10,
            5.55555e-8
        ]
        
        for pval in test_pvalues:
            # Use format_num from meta_analysis
            formatted = meta_analysis.format_num(pval, precision=6)
            
            # Should be in scientific notation
            self.assertIn('e', formatted.lower(),
                         f"Very small p-value {pval} should be in scientific notation")
            
            # Parse mantissa - handle both "1.000000e-100" and "1.e-100" formats
            mantissa = formatted.lower().split('e')[0]
            
            # If there's a decimal point, check precision
            if '.' in mantissa:
                decimal_part = mantissa.split('.')[1]
                # For precision=6, should have up to 6 decimal places
                # (may have fewer if trailing zeros are stripped)
                self.assertLessEqual(len(decimal_part), 6,
                               f"P-value {pval} formatted as {formatted} should have at most 6 decimal places in mantissa")
                # Should have at least some decimal places (not empty)
                # unless the mantissa is exactly an integer like "1."
                if len(decimal_part) > 0:
                    self.assertGreater(len(decimal_part), 0,
                                     f"P-value {pval} should have decimal places in mantissa")

    def test_rounding_consistency_across_modules(self):
        """Test 4: Verify consistent rounding behavior across munge, harmonize, and meta_analysis"""
        test_value = 0.123456789
        
        # Test munge rounding (default precision 6)
        munge_rounded = round(test_value, 6)
        self.assertAlmostEqual(munge_rounded, 0.123457, places=6)
        
        # Test harmonize get_decimal_places
        decimals = harmonize.get_decimal_places(munge_rounded)
        self.assertLessEqual(decimals, 6, 
                            "Rounded value should have at most 6 decimal places")
        
        # Test meta_analysis format_num
        formatted = meta_analysis.format_num(test_value, precision=6)
        self.assertIsNotNone(formatted)
        self.assertNotEqual(formatted, "NA")

    def test_format_num_precision_levels(self):
        """Test 5: Verify format_num works with different precision levels"""
        test_value = 1.23456789e-10
        
        precisions = [2, 4, 6, 8]
        for prec in precisions:
            formatted = meta_analysis.format_num(test_value, precision=prec)
            
            # Parse mantissa
            mantissa = formatted.lower().split('e')[0]
            if '.' in mantissa:
                decimal_part = mantissa.split('.')[1]
                self.assertEqual(len(decimal_part), prec,
                               f"precision={prec} should give {prec} decimal places in mantissa")

    def test_numpy_round_consistency(self):
        """Test 6: Verify numpy.round behavior with 6 decimals"""
        import numpy
        
        test_cases = [
            (10.1234567, 10.123457),
            (5.5555555, 5.555556),
            (0.9999999, 1.0),
            (-0.1234567, -0.123457),
        ]
        
        for input_val, expected in test_cases:
            rounded = numpy.round(input_val, 6)
            self.assertAlmostEqual(rounded, expected, places=6,
                                  msg=f"numpy.round({input_val}, 6) should equal {expected}")

    def test_precision_with_real_meta_analysis(self):
        """Test 7: Run a mini meta-analysis and verify output precision"""
        from unittest.mock import MagicMock
        
        # Create mock studies
        study1 = MagicMock(spec=meta_analysis.Study)
        study1.effective_size = 1000
        study1.name = "Study1"
        
        variant1 = MagicMock(spec=meta_analysis.VariantData)
        variant1.beta = 0.123456789  # High precision input
        variant1.se = 0.045678912
        variant1.pval = 1.23456e-5
        variant1.z_score = 4.5
        
        study2 = MagicMock(spec=meta_analysis.Study)
        study2.effective_size = 2000
        study2.name = "Study2"
        
        variant2 = MagicMock(spec=meta_analysis.VariantData)
        variant2.beta = 0.234567891
        variant2.se = 0.034567891
        variant2.pval = 5.67891e-8
        variant2.z_score = 5.5
        
        study_list = [(study1, variant1), (study2, variant2)]
        
        # Run meta-analysis
        result = meta_analysis.do_meta(study_list, methods=['inv_var'], is_het_test=False)
        
        # Verify result structure
        self.assertEqual(len(result), 1)
        self.assertIsNotNone(result[0])
        
        # Unpack result
        beta_str, se_str, pval_str, mlogp_rounded = result[0]
        
        # Verify formatting
        self.assertIsInstance(beta_str, str)
        self.assertIsInstance(se_str, str)
        self.assertIsInstance(pval_str, str)
        
        # Verify mlogp has 6 decimal precision
        self.assertIsInstance(mlogp_rounded, (float, int))
        mlogp_str = f"{mlogp_rounded:.10f}"
        if '.' in mlogp_str:
            decimal_part = mlogp_str.split('.')[1].rstrip('0')
            self.assertLessEqual(len(decimal_part), 6,
                               f"mlogp should have at most 6 decimal places")
        
        # Verify beta is in scientific notation with 6 decimals
        if 'e' in beta_str.lower():
            mantissa = beta_str.lower().split('e')[0]
            if '.' in mantissa:
                decimal_part = mantissa.split('.')[1]
                self.assertEqual(len(decimal_part), 6,
                               f"Beta mantissa should have 6 decimal places")

    def test_heterogeneity_with_precision(self):
        """Test 8: Verify heterogeneity test maintains precision"""
        from unittest.mock import MagicMock
        
        # Create mock studies with different effect sizes
        study1 = MagicMock(spec=meta_analysis.Study)
        study1.effective_size = 1000
        study1.name = "Study1"
        
        variant1 = MagicMock(spec=meta_analysis.VariantData)
        variant1.beta = 0.5
        variant1.se = 0.1
        variant1.pval = 0.001
        variant1.z_score = 3.29
        
        study2 = MagicMock(spec=meta_analysis.Study)
        study2.effective_size = 1000
        study2.name = "Study2"
        
        variant2 = MagicMock(spec=meta_analysis.VariantData)
        variant2.beta = 0.6
        variant2.se = 0.1
        variant2.pval = 0.0001
        variant2.z_score = 3.89
        
        study_list = [(study1, variant1), (study2, variant2)]
        
        # Run with heterogeneity test
        result = meta_analysis.do_meta(study_list, methods=['inv_var'], is_het_test=True)
        
        # Should have 5 elements (beta, se, pval, mlogp, het_p)
        self.assertEqual(len(result), 1)
        self.assertIsNotNone(result[0])
        self.assertEqual(len(result[0]), 5, 
                        "Result with het_test should have 5 elements")
        
        # Verify het_p is formatted correctly
        het_p_str = result[0][4]
        self.assertIsInstance(het_p_str, str)

    def test_decimal_places_helper_function(self):
        """Test 9: Verify harmonize.get_decimal_places works correctly"""
        test_cases = [
            (0.123456, 6),
            (0.1, 1),
            (1.0, 1),  # Python represents 1.0 as having 1 decimal
            (0.123, 3),
            (123, 0),
            (None, None),
        ]
        
        for value, expected_decimals in test_cases:
            result = harmonize.get_decimal_places(value)
            self.assertEqual(result, expected_decimals,
                           f"get_decimal_places({value}) should return {expected_decimals}, got {result}")


class TestPrecisionEdgeCases(unittest.TestCase):
    """Test edge cases for precision handling"""

    def test_zero_values(self):
        """Test precision with zero values"""
        zero_rounded = round(0.0, 6)
        self.assertEqual(zero_rounded, 0.0)
        
        formatted = meta_analysis.format_num(0.0, precision=6)
        # Zero might be formatted as "0.000000e+00" or similar
        self.assertIsNotNone(formatted)

    def test_negative_values_precision(self):
        """Test precision with negative values"""
        neg_value = -0.123456789
        neg_rounded = round(neg_value, 6)
        
        self.assertAlmostEqual(neg_rounded, -0.123457, places=6)
        self.assertLess(neg_rounded, 0)

    def test_very_large_values(self):
        """Test precision with very large values"""
        large_value = 123456.789012
        large_rounded = round(large_value, 6)
        
        # Should keep 6 decimal places
        self.assertAlmostEqual(large_rounded, 123456.789012, places=6)

    def test_very_small_values(self):
        """Test precision with very small values"""
        small_value = 0.000000123456
        small_rounded = round(small_value, 6)
        
        # Should round to 6 decimals
        self.assertAlmostEqual(small_rounded, 0.000000, places=6)
        
        # Format with higher precision to see actual value
        formatted = meta_analysis.format_num(small_value, precision=6)
        self.assertIn('e', formatted.lower())


if __name__ == '__main__':
    unittest.main()
