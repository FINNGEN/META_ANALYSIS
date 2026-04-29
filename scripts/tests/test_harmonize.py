#!/usr/bin/env python3

import unittest
import sys
import os
import tempfile
import gzip
from unittest.mock import patch, MagicMock

# Add parent directory to path to import harmonize
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harmonize


class TestHarmonizePrecision(unittest.TestCase):
    """Test suite for validating numeric precision preservation in harmonize.py"""

    def test_get_decimal_places_6_decimals(self):
        """Test 1: Verify get_decimal_places correctly counts 6 decimals"""
        result = harmonize.get_decimal_places(0.123456)
        self.assertEqual(result, 6, 
                        "get_decimal_places should return 6 for a number with 6 decimal places")

    def test_get_decimal_places_more_than_6(self):
        """Test 2: Verify get_decimal_places handles >6 decimals"""
        result = harmonize.get_decimal_places(0.123456789)
        self.assertEqual(result, 9, 
                        "get_decimal_places should return 9 for a number with 9 decimal places")

    def test_get_decimal_places_integer(self):
        """Test 3: Verify get_decimal_places handles integers"""
        result = harmonize.get_decimal_places(123)
        self.assertEqual(result, 0, 
                        "get_decimal_places should return 0 for an integer")

    def test_get_decimal_places_none(self):
        """Test 4: Verify get_decimal_places handles None"""
        result = harmonize.get_decimal_places(None)
        self.assertIsNone(result, 
                         "get_decimal_places should return None when input is None")

    def test_flip_beta_preserves_precision(self):
        """Test 5: Verify beta negation maintains precision"""
        beta_original = 0.123456
        decimal_places = harmonize.get_decimal_places(beta_original)
        
        # Simulate flipping
        beta_flipped = -1 * beta_original
        beta_flipped_rounded = round(beta_flipped, decimal_places)
        
        # Check precision is maintained
        self.assertEqual(harmonize.get_decimal_places(beta_flipped_rounded), decimal_places,
                        "Flipped beta should maintain the same number of decimal places")
        self.assertAlmostEqual(beta_flipped_rounded, -0.123456, places=6,
                              msg="Flipped beta should preserve precision")

    def test_flip_af_preserves_precision(self):
        """Test 6: Verify AF complement (1-af) maintains precision"""
        af_original = 0.234567
        decimal_places = harmonize.get_decimal_places(af_original)
        
        # Simulate AF flipping (1 - af)
        af_flipped = 1 - af_original
        af_flipped_rounded = round(af_flipped, decimal_places)
        
        # Check precision is maintained
        result_decimals = harmonize.get_decimal_places(af_flipped_rounded)
        self.assertLessEqual(result_decimals, decimal_places + 1,  # Allow for minor floating point
                            "Flipped AF should maintain similar decimal places")
        self.assertAlmostEqual(af_flipped_rounded, 0.765433, places=6,
                              msg="Flipped AF should preserve precision")

    def test_precision_with_scientific_notation(self):
        """Test 7: Verify precision handling with very small numbers"""
        # Very small number that might be in scientific notation
        small_number = 1.23456e-8
        
        # Convert to string and count decimal places
        str_num = f"{small_number:.10f}"  # Format with enough precision
        if '.' in str_num:
            # Count trailing non-zero decimals
            decimal_part = str_num.split('.')[1].rstrip('0')
            decimal_places = len(decimal_part)
            self.assertGreaterEqual(decimal_places, 6,
                                   "Small numbers should maintain at least 6 significant digits")

    def test_rounding_consistency(self):
        """Test 8: Verify consistent rounding behavior"""
        # Test that rounding to detected decimal places is consistent
        numbers = [0.123456, 0.987654, 0.555555, 0.111111]
        
        for num in numbers:
            places = harmonize.get_decimal_places(num)
            rounded = round(num, places)
            
            # Should round to itself if already at that precision
            self.assertAlmostEqual(rounded, num, places=places,
                                  msg=f"Number {num} should round to itself at {places} decimal places")

    def test_decimal_places_with_trailing_zeros(self):
        """Test 9: Verify handling of numbers with trailing zeros"""
        # 0.100000 as a string has 6 decimals, but as float it's 0.1
        num = 0.100000
        places = harmonize.get_decimal_places(num)
        
        # Python will strip trailing zeros, so this should be 1
        self.assertEqual(places, 1,
                        "get_decimal_places should not count trailing zeros in float representation")

    def test_large_decimal_precision(self):
        """Test 10: Verify handling of very high precision numbers"""
        # Number with 12 decimal places
        high_precision = 0.123456789012
        places = harmonize.get_decimal_places(high_precision)
        
        self.assertEqual(places, 12,
                        "get_decimal_places should correctly count 12 decimal places")
        
        # Round to 6 and verify
        rounded_6 = round(high_precision, 6)
        places_6 = harmonize.get_decimal_places(rounded_6)
        self.assertLessEqual(places_6, 6,
                            "After rounding to 6, should have at most 6 decimal places")


class TestHarmonizeFunctionality(unittest.TestCase):
    """Test harmonize function behavior (focusing on precision)"""

    def setUp(self):
        """Set up test fixtures"""
        # Create temporary test files
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up test files"""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_harmonize_preserves_6_decimal_precision(self):
        """Test: Verify harmonization preserves 6 decimal place values"""
        # This is a placeholder test - actual implementation would require
        # creating mock input files and running harmonize function
        
        # Test the principle: if we have 6 decimal places, rounding should preserve them
        beta_6_decimals = 0.123456
        places = harmonize.get_decimal_places(beta_6_decimals)
        rounded = round(beta_6_decimals, places)
        
        self.assertEqual(harmonize.get_decimal_places(rounded), 6,
                        "6 decimal place values should be preserved through rounding")


if __name__ == '__main__':
    unittest.main()
