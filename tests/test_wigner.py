#!/usr/bin/env python3
"""
Tests for the Wigner class in the Wigner6j package.

The Wigner class represents a complete Wigner 6j-symbol with six Young diagrams
and four vertices, providing methods to validate and compute the symbol.
"""
import unittest
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from wigner_new import Wigner
from wigner_new import Fraction

class TestWigner(unittest.TestCase):
    """Test cases for the Wigner class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Common Young diagrams for testing
        self.trivial = [0]
        self.quark = [1, 0, 0]
        self.antiquark = [0, 0, 1]
        self.gluon = [2, 1, 0]  # Adjoint representation for SU(3)
        self.complex_diagram = [3, 2, 1]
    
    def test_initialise_simple(self):
        """Test Wigner Initialisation with simple diagrams."""

        #Case 0
        self.assertEqual(Wigner((2),(3),(1),(2),(1),(1),1,1,1,1,3).get_value(),Fraction(1,6,1,1))
        self.assertEqual(Wigner((3,1),(2),(1),(1),(2,1),(1),1,1,1,1,3).get_value(),Fraction(1,1,1,48).reduce())
        self.assertEqual(Wigner((1),(1,1),(1),(1),(2,1),(1),1,1,1,1,3).get_value(),Fraction(1,6,1,1))
        
        
        
    def test_case_1_symbols(self):# Case 1
        """Test Wigner 6j-symbols of Case 1."""
        
        # Vertex 1 alpha=gamma and real, vertex 2 is 1+
        self.assertEqual(Wigner((2,1),(2,2),(2,1),(1),(2,1),(1),1,1,1,1,3).get_value(),Fraction(-1,8,5,6).reduce())
        self.assertEqual(Wigner((4,2),(4,3),(4,2),(1),(2,1),(1),1,1,1,1,3).get_value(),Fraction(-1,72,35,2).reduce())
        self.assertEqual(Wigner((4,2),(5,2),(4,2),(1),(2,1),(1),1,1,1,1,3).get_value(),Fraction(1,36,5,14).reduce())
        
        # Vertex 1 alpha=gamma and real, vertex 2 is 1-
        self.assertEqual(Wigner((2,1),(2,2),(2,1),(1),(2,1),(1),1,2,1,1,3).get_value(),Fraction(-1,8,1,6).reduce())
        self.assertEqual(Wigner((4,2),(4,3),(4,2),(1),(2,1),(1),1,2,1,1,3).get_value(),Fraction(-1,24,1,6).reduce())
        self.assertEqual(Wigner((4,2),(5,2),(4,2),(1),(2,1),(1),1,2,1,1,3).get_value(),Fraction(1,12,1,6).reduce())
    
    def test_case_1_symbols_2(self):# Case 1
        """Test Wigner 6j-symbols of Case 1."""
        
        # Vertex 1 alpha=gamma and not real
        self.assertEqual(Wigner((3,1),(3,2),(3,1),(1),(2,1),(1),1,1,1,1,3).get_value(),Fraction(-1,3,1,35).reduce())
        self.assertEqual(Wigner((3,1),(4,1),(3,1),(1),(2,1),(1),1,1,1,1,3).get_value(),Fraction(1,24,7,5).reduce())
        self.assertEqual(Wigner((4,1),(5,1),(4,1),(1),(2,1),(1),1,1,1,1,3).get_value(),Fraction(1,24,37,35).reduce())

        self.assertEqual(Wigner((3,1),(3,2),(3,1),(1),(2,1),(1),1,2,1,1,3).get_value(),Fraction(1,2,1,105).reduce())
        self.assertEqual(Wigner((3,1),(4,1),(3,1),(1),(2,1),(1),1,2,1,1,3).get_value(),Fraction(0).reduce())
        self.assertEqual(Wigner((4,1),(5,1),(4,1),(1),(2,1),(1),1,2,1,1,3).get_value(),Fraction(0).reduce())
        
        self.assertEqual(Wigner((4,1),(4,2),(4,1),(1),(2,1),(1),1,2,1,1,3).get_value(),Fraction(1,6,5,111).reduce())



if __name__ == '__main__':
    unittest.main()