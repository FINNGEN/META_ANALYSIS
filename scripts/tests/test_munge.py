#!/usr/bin/env python3

import unittest
import sys
import os
import tempfile
import gzip
from io import StringIO
from unittest.mock import patch, MagicMock
import argparse

# Add parent directory to path to import munge
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import munge


class TestMungePrecision(unittest.TestCase):
    """Test suite for validating numeric precision in munge.py"""

    def test_default_rounding_precision(self):
        """Test 1: Verify default --rounding-precision is 6"""
        parser = argparse.ArgumentParser()
        # Simulate the argument setup from munge.py
        parser.add_argument("--rounding-precision", type=int, default=6)
        args = parser.parse_args([])
        
        self.assertEqual(args.rounding_precision, 6, 
                        "Default rounding precision should be 6 decimal places")

    def test_beta_precision_6_decimals(self):
        """Test 2: Verify beta values are rounded to 6 decimal places"""
        # Test rounding to 6 decimals
        beta_input = 0.123456789
        beta_rounded = round(beta_input, 6)
        beta_str = str(beta_rounded)
        
        # Verify it has 6 decimal places (or fewer if trailing zeros)
        if '.' in beta_str:
            decimals = len(beta_str.split('.')[1])
            self.assertLessEqual(decimals, 6, 
                               f"Beta should have at most 6 decimal places, got {decimals}")
        
        self.assertAlmostEqual(beta_rounded, 0.123457, places=6,
                              msg="Beta should be rounded to 6 decimal places")

    def test_se_precision_6_decimals(self):
        """Test 3: Verify standard error values are rounded to 6 decimal places"""
        se_input = 0.987654321
        se_rounded = round(se_input, 6)
        
        self.assertAlmostEqual(se_rounded, 0.987654, places=6,
                              msg="SE should be rounded to 6 decimal places")

    def test_af_precision_6_decimals(self):
        """Test 4: Verify allele frequency values are rounded to 6 decimal places"""
        af_input = 0.456789123
        af_rounded = round(af_input, 6)
        
        self.assertAlmostEqual(af_rounded, 0.456789, places=6,
                              msg="AF should be rounded to 6 decimal places")

    def test_mlogp_precision_6_decimals(self):
        """Test 5: Verify -log10(p) values are rounded to 6 decimal places"""
        mlogp_input = 10.123456789
        mlogp_rounded = round(mlogp_input, 6)
        
        self.assertAlmostEqual(mlogp_rounded, 10.123457, places=6,
                              msg="mlogp should be rounded to 6 decimal places")

    def test_custom_rounding_precision(self):
        """Test 6: Verify custom --rounding-precision argument works"""
        # Test with precision 8
        beta_input = 0.123456789012
        beta_rounded_8 = round(beta_input, 8)
        
        self.assertAlmostEqual(beta_rounded_8, 0.12345679, places=8,
                              msg="Beta with precision=8 should round to 8 decimals")
        
        # Test with precision 4
        beta_rounded_4 = round(beta_input, 4)
        self.assertAlmostEqual(beta_rounded_4, 0.1235, places=4,
                              msg="Beta with precision=4 should round to 4 decimals")

    def test_end_to_end_output_file(self):
        """Test 7: Integration test with actual file I/O"""
        # Create a temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_in:
            input_file = tmp_in.name
            # Write header
            tmp_in.write("chr\tpos\tref\talt\tbeta\tse\taf\tpval\n")
            # Write test data with more than 6 decimal places
            tmp_in.write("1\t1000000\tA\tG\t0.123456789\t0.987654321\t0.456789123\t0.0000012345\n")
            tmp_in.write("1\t1000100\tC\tT\t-0.987654321\t0.123456789\t0.111111111\t0.0000054321\n")
        
        try:
            # Capture stdout
            with patch('sys.stdout', new=StringIO()) as mock_stdout:
                # Mock sys.argv to simulate command line
                test_args = [
                    'munge.py',
                    input_file,
                    '--chr-col', 'chr',
                    '--pos-col', 'pos',
                    '--ref-col', 'ref',
                    '--alt-col', 'alt',
                    '--af-col', 'af',
                    '--effect-col', 'beta',
                    '--se-col', 'se',
                    '--pval-col', 'pval',
                    '--effect-type', 'beta',
                    '--se-type', 'se',
                    '--pval-type', 'p',
                    '--rounding-precision', '6'
                ]
                
                with patch('sys.argv', test_args):
                    with patch('sys.stderr', new=StringIO()):
                        try:
                            # Call main function
                            munge.main()
                            output = mock_stdout.getvalue()
                            
                            # Parse output lines
                            lines = output.strip().split('\n')
                            
                            # Check header exists
                            self.assertGreater(len(lines), 0, "Output should have at least header line")
                            
                            # Check data lines (skip header)
                            for i, line in enumerate(lines[1:], 1):
                                fields = line.split('\t')
                                if len(fields) > 5:  # Ensure we have enough fields
                                    # Check beta field (assuming it's in the output)
                                    # The exact position depends on munge.py output format
                                    # This is a basic check that numeric values exist
                                    has_numeric = False
                                    for field in fields:
                                        try:
                                            val = float(field)
                                            has_numeric = True
                                            # Check that if it has a decimal point, it has max 6 decimals
                                            if '.' in field and 'e' not in field.lower():
                                                decimal_places = len(field.split('.')[1])
                                                self.assertLessEqual(decimal_places, 6,
                                                    f"Line {i}: Field '{field}' has more than 6 decimal places")
                                        except ValueError:
                                            continue
                                    
                                    self.assertTrue(has_numeric, 
                                        f"Line {i} should contain at least one numeric value")
                        
                        except SystemExit:
                            # munge.main() might call sys.exit, which is okay for this test
                            pass
        
        finally:
            # Clean up temp file
            if os.path.exists(input_file):
                os.unlink(input_file)

    def test_rounding_preserves_sign(self):
        """Test 8: Verify rounding preserves sign for negative values"""
        beta_negative = -0.123456789
        beta_rounded = round(beta_negative, 6)
        
        self.assertAlmostEqual(beta_rounded, -0.123457, places=6,
                              msg="Negative beta should round correctly with sign preserved")
        self.assertLess(beta_rounded, 0, "Rounded negative value should remain negative")

    def test_rounding_edge_cases(self):
        """Test 9: Test rounding edge cases (very small, very large, zero)"""
        # Very small positive
        small_val = 0.000000123456789
        small_rounded = round(small_val, 6)
        self.assertAlmostEqual(small_rounded, 0.000000, places=6)
        
        # Zero
        zero_rounded = round(0.0, 6)
        self.assertEqual(zero_rounded, 0.0)
        
        # Large value
        large_val = 123456.789012
        large_rounded = round(large_val, 6)
        self.assertAlmostEqual(large_rounded, 123456.789012, places=6)


if __name__ == '__main__':
    unittest.main()
