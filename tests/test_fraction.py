#!/usr/bin/env python3
"""
Tests for the Fraction class in the Wigner6j package.

The Fraction class handles exact arithmetic without floating point errors,
supporting fractions with square root terms in numerator and denominator.
"""
import unittest
import math
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from wigner_new import Fraction


class TestFraction(unittest.TestCase):
    """Test cases for the Fraction class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.zero = Fraction(0, 1)
        self.one = Fraction(1, 1)
        self.half = Fraction(1, 2)
        self.two_thirds = Fraction(2, 3)
        self.three_halves = Fraction(3, 2)
        self.negative_half = Fraction(-1, 2)
    
    def test_initialization(self):
        """Test Fraction initialization with different parameters."""
        # Basic fraction
        f = Fraction(3, 4)
        self.assertEqual(f.nominator, 3)
        self.assertEqual(f.denominator, 4)
        self.assertEqual(f.nominatorRoot, 1)
        self.assertEqual(f.denominatorRoot, 1)
        
        # Fraction with roots
        f = Fraction(2, 3, 5, 7)
        self.assertEqual(f.nominator, 2)
        self.assertEqual(f.denominator, 3)
        self.assertEqual(f.nominatorRoot, 5)
        self.assertEqual(f.denominatorRoot, 7)
        
        # Default values
        f = Fraction()
        self.assertEqual(f.nominator, 0)
        self.assertEqual(f.denominator, 1)
        self.assertEqual(f.nominatorRoot, 1)
        self.assertEqual(f.denominatorRoot, 1)
    
    def test_initialization_type_conversion(self):
        """Test that non-integer inputs are properly converted."""
        # Float inputs should be rounded to integers
        f = Fraction(3.7, 4.2, 5.1, 7.9)
        self.assertEqual(f.nominator, 4)
        self.assertEqual(f.denominator, 4)
        self.assertEqual(f.nominatorRoot, 5)
        self.assertEqual(f.denominatorRoot, 8)
    
    def test_setattr_type_checking(self):
        """Test that setting attributes enforces integer types."""
        f = Fraction(1, 2)
        with self.assertRaises(TypeError):
            f.nominator = 3.5
        with self.assertRaises(TypeError):
            f.denominator = "invalid"
    
    def test_addition(self):
        """Test fraction addition."""
        # Basic addition
        result = self.half + self.two_thirds
        # 1/2 + 2/3 = 3/6 + 4/6 = 7/6
        self.assertEqual(result.nominator, 7)
        self.assertEqual(result.denominator, 6)
        
        # Addition with zero
        result = self.half + self.zero
        self.assertEqual(result.nominator, 1)
        self.assertEqual(result.denominator, 2)
        
        # Addition with integer
        result = self.half + 2
        # 1/2 + 2 = 1/2 + 4/2 = 5/2
        self.assertEqual(result.nominator, 5)
        self.assertEqual(result.denominator, 2)
        
        # Right addition with integer
        result = 2 + self.half
        self.assertEqual(result.nominator, 5)
        self.assertEqual(result.denominator, 2)
    
    def test_subtraction(self):
        """Test fraction subtraction."""
        # Basic subtraction
        result = self.two_thirds - self.half
        # 2/3 - 1/2 = 4/6 - 3/6 = 1/6
        self.assertEqual(result.nominator, 1)
        self.assertEqual(result.denominator, 6)
        
        # Subtraction with integer
        result = self.three_halves - 1
        # 3/2 - 1 = 3/2 - 2/2 = 1/2
        self.assertEqual(result.nominator, 1)
        self.assertEqual(result.denominator, 2)
    
    def test_multiplication(self):
        """Test fraction multiplication."""
        # Basic multiplication
        result = self.half * self.two_thirds
        # 1/2 * 2/3 = 2/6 = 1/3
        result = result.reduce()
        self.assertEqual(result.nominator, 1)
        self.assertEqual(result.denominator, 3)
        
        # Multiplication with integer
        result = self.half * 3
        self.assertEqual(result.nominator, 3)
        self.assertEqual(result.denominator, 2)
        
        # Right multiplication with integer
        result = 3 * self.half
        self.assertEqual(result.nominator, 3)
        self.assertEqual(result.denominator, 2)
    
    def test_division(self):
        """Test fraction division."""
        # Basic division
        result = self.half / self.two_thirds
        # 1/2 ÷ 2/3 = 1/2 * 3/2 = 3/4
        result = result.reduce()
        self.assertEqual(result.nominator, 3)
        self.assertEqual(result.denominator, 4)
        
        # Division by integer
        result = self.half / 2
        # 1/2 ÷ 2 = 1/4
        self.assertEqual(result.nominator, 1)
        self.assertEqual(result.denominator, 4)
        
        # Integer divided by fraction
        result = 3 / self.half
        # 3 ÷ 1/2 = 3 * 2/1 = 6
        self.assertEqual(result.nominator, 6)
        self.assertEqual(result.denominator, 1)
    
    def test_comparison_operators(self):
        """Test comparison operators."""
        # Equality
        self.assertTrue(self.half == Fraction(1, 2))
        self.assertTrue(self.half == Fraction(2, 4))  # After reduction
        self.assertFalse(self.half == self.two_thirds)
        
        # Less than
        self.assertTrue(self.half < self.two_thirds)
        self.assertFalse(self.two_thirds < self.half)
        
        # Greater than
        self.assertTrue(self.two_thirds > self.half)
        self.assertFalse(self.half > self.two_thirds)
        
        # Comparison with integers/floats
        self.assertTrue(self.half < 1)
        self.assertTrue(self.half > 0)
        self.assertTrue(self.half == 0.5)
    
    def test_reduce(self):
        """Test fraction reduction/simplification."""
        # Test basic reduction
        f = Fraction(6, 9)
        f = f.reduce()
        self.assertEqual(f.nominator, 2)
        self.assertEqual(f.denominator, 3)
        
        # Test negative denominator handling
        f = Fraction(1, -2)
        f = f.reduce()
        self.assertEqual(f.nominator, -1)
        self.assertEqual(f.denominator, 2)
        
        # Test zero denominator
        f = Fraction(1, 0)
        result = f.reduce()
        self.assertIsNone(result)
    
    def test_sqrt(self):
        """Test square root operation."""
        f = Fraction(4, 9)
        sqrt_f = f.sqrt()
        # sqrt(4/9) should give us a fraction with roots
        self.assertEqual(sqrt_f.nominatorRoot, 4)
        self.assertEqual(sqrt_f.denominatorRoot, 9)
    
    def test_square(self):
        """Test squaring operation."""
        f = Fraction(2, 3)
        squared = f.sq()
        self.assertEqual(squared.nominator, 4)
        self.assertEqual(squared.denominator, 9)
    
    def test_power(self):
        """Test power operation."""
        f = Fraction(2, 3)
        cubed = f ** 3
        self.assertEqual(cubed.nominator, 8)
        self.assertEqual(cubed.denominator, 27)
    
    def test_absolute_value(self):
        """Test absolute value operation."""
        f = Fraction(-3, 4)
        abs_f = f.abs()
        self.assertEqual(abs_f.nominator, 3)
        self.assertEqual(abs_f.denominator, 4)
        
        # Test with positive value
        f = Fraction(3, 4)
        abs_f = f.abs()
        self.assertEqual(abs_f.nominator, 3)
        self.assertEqual(abs_f.denominator, 4)
    
    def test_string_representation(self):
        """Test string representation methods."""
        # Test __str__
        self.assertEqual(str(self.zero), "0")
        self.assertEqual(str(self.one), "1")
        self.assertEqual(str(self.half), "1/2")
        
        # Test __repr__
        self.assertEqual(repr(self.zero), "0")
        self.assertEqual(repr(self.half), "1/2")
    
    def test_edge_cases(self):
        """Test edge cases and error conditions."""
        # Division by zero should be handled gracefully
        self.assertEqual(self.one / self.zero, None)
    
    def test_root_operations(self):
        """Test operations involving square roots."""
        # Create fractions with roots
        f1 = Fraction(1, 1, 2, 1)  # (1*sqrt(2))/1
        f2 = Fraction(2, 1, 1, 3)  # 2/sqrt(3)
        
        # Test that operations preserve root structure
        result = f1 * f2
        # Should have roots in both numerator and denominator
        self.assertTrue(result.nominatorRoot != 1 or result.denominatorRoot != 1)
    
    def test_precision_preservation(self):
        """Test that exact arithmetic is preserved."""
        # Operations that would lose precision in floating point
        # should be exact with fractions
        f = Fraction(1, 3)
        result = f * 3
        result = result.reduce()
        self.assertEqual(result.nominator, 1)
        self.assertEqual(result.denominator, 1)  # Should be exactly 1
    
    def test_zero_handling(self):
        """Test special behavior with zero values."""
        # Zero in numerator
        zero_frac = Fraction(0, 5)
        self.assertEqual(str(zero_frac), "0")
        
        # Operations with zero
        result = self.half + zero_frac
        self.assertEqual(result.nominator, 1)
        self.assertEqual(result.denominator, 2)
        
        result = self.half * zero_frac
        self.assertEqual(result.nominator, 0)


if __name__ == '__main__':
    unittest.main()