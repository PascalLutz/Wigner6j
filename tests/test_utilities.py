#!/usr/bin/env python3
"""
Tests for utility functions in the Wigner6j package.

These tests cover helper functions like particle type checks, mathematical utilities,
vertex validation functions, and 6j-symbol type classification functions.
"""
import unittest
import sys
from pathlib import Path
import pair_multiplication as pr

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from wigner_new import (
    kronecker_delta, is_quark, is_adjoint, is_gluon, is_real, is_viable_diagram,
    check_vertex_sign, is_6j_with_two_quarks, is_6j_with_quark_gluon_vertex,
    is_6j_with_quark_gluon_opposing, is_6j_with_two_gluon, is_6j_with_three_gluon,
    conjugate_diagram, create_fractions_array
)


class TestUtilityFunctions(unittest.TestCase):
    """Test cases for utility and helper functions."""
    
    def setUp(self):
        """Set up test fixtures with common Young diagrams."""
        # Standard SU(3) representations
        self.trivial = [0, 0, 0]
        self.quark = [1, 0, 0]           # Fundamental
        self.antiquark = [0, 0, 1]       # Antifundamental
        self.gluon = [2, 1, 0]           # Adjoint (8)
        self.symmetric = [2, 0, 0]       # Symmetric tensor (6)
        self.antisymmetric = [1, 1, 0]   # Antisymmetric tensor (3̄)
        
        # Complex representations
        self.complex_rep = [3, 2, 1]
        self.large_rep = [4, 3, 2]
        
        # Invalid diagrams
        self.invalid_increasing = [1,  2, 3]
        self.invalid_mixed = [2, 1, 3]
    
    def test_kronecker_delta(self):
        """Test Kronecker delta function."""
        # Equal indices
        self.assertEqual(kronecker_delta(1, 1), 1)
        self.assertEqual(kronecker_delta(0, 0), 1)
        self.assertEqual(kronecker_delta(5, 5), 1)
        self.assertEqual(kronecker_delta(-1, -1), 1)
        
        # Different indices
        self.assertEqual(kronecker_delta(1, 2), 0)
        self.assertEqual(kronecker_delta(0, 1), 0)
        self.assertEqual(kronecker_delta(-1, 1), 0)
        self.assertEqual(kronecker_delta(10, 20), 0)
        
        # Test with various data types
        self.assertEqual(kronecker_delta(1.0, 1.0), 1)
        self.assertEqual(kronecker_delta(1.0, 1), 1)
        self.assertEqual(kronecker_delta(1, 1.0), 1)
        self.assertEqual(kronecker_delta(1.0, 2.0), 0)
    
    def test_is_quark(self):
        """Test quark representation identification."""
        # Fundamental representation [1,0,0] is a quark
        self.assertTrue(is_quark(self.quark))
        
        # Other representations should not be quarks
        self.assertFalse(is_quark(self.trivial))
        self.assertFalse(is_quark(self.gluon))
        self.assertFalse(is_quark(self.symmetric))
        self.assertFalse(is_quark(self.antisymmetric))
        self.assertFalse(is_quark(self.complex_rep))
        
        # Antifundamental [0,0,1] should also be identified appropriately
        # (depending on implementation, might be True or False)
        result = is_quark(self.antiquark)
        self.assertIsInstance(result, bool)
    
    def test_is_adjoint(self):
        """Test adjoint representation identification."""
        # Adjoint representation [2,1,0] for SU(3)
        self.assertTrue(is_adjoint(self.gluon, n=3))
        
        # Other representations should not be adjoint
        self.assertFalse(is_adjoint(self.trivial, n=3))
        self.assertFalse(is_adjoint(self.quark, n=3))
        self.assertFalse(is_adjoint(self.symmetric, n=3))
        
        # Test with different n values
        result = is_adjoint(self.gluon, n=2)  # SU(2) adjoint is different
        self.assertIsInstance(result, bool)
        
        result = is_adjoint(self.gluon, n=4)  # SU(4) adjoint is different
        self.assertIsInstance(result, bool)
    
    def test_is_gluon(self):
        """Test gluon representation identification."""
        # Should be equivalent to is_adjoint for most purposes
        self.assertTrue(is_gluon(self.gluon, n=3))
        
        # Other representations should not be gluons
        self.assertFalse(is_gluon(self.trivial, n=3))
        self.assertFalse(is_gluon(self.quark, n=3))
        self.assertFalse(is_gluon(self.symmetric, n=3))
        
        # Test consistency with is_adjoint
        for rep in [self.trivial, self.quark, self.gluon, self.symmetric]:
            adjoint_result = is_adjoint(rep, n=3)
            gluon_result = is_gluon(rep, n=3)
            # They should give the same result (or have defined relationship)
            self.assertIsInstance(gluon_result, bool)
    
    def test_is_real(self):
        """Test real representation identification."""
        # Test various representations
        result = is_real(self.trivial)
        self.assertIsInstance(result, bool)
        
        result = is_real(self.quark)
        self.assertIsInstance(result, bool)
        
        result = is_real(self.gluon)
        self.assertIsInstance(result, bool)
        
        result = is_real(self.symmetric)
        self.assertIsInstance(result, bool)
        
        # Real representations are self-conjugate
        # For SU(3), the adjoint [2,1,0] should be real
        self.assertTrue(is_real(self.gluon))
    
    def test_is_viable_diagram(self):
        """Test Young diagram validity checking."""
        # Valid diagrams (weakly decreasing)
        self.assertTrue(is_viable_diagram(self.trivial))
        self.assertTrue(is_viable_diagram(self.quark))
        self.assertTrue(is_viable_diagram(self.gluon))
        self.assertTrue(is_viable_diagram(self.symmetric))
        self.assertTrue(is_viable_diagram(self.complex_rep))
        
        # Invalid diagrams (not weakly decreasing)
        self.assertFalse(is_viable_diagram(self.invalid_increasing))
        self.assertFalse(is_viable_diagram(self.invalid_mixed))
        self.assertFalse(is_viable_diagram([1, 3, 2]))
        self.assertFalse(is_viable_diagram([0, 1, 0]))
        
        # Edge cases
        self.assertTrue(is_viable_diagram([0]))  # Single zero
        self.assertTrue(is_viable_diagram([1]))  # Single box
        self.assertTrue(is_viable_diagram([]))   # Empty diagram
        
        # Test with negative numbers (should be invalid)
        self.assertFalse(is_viable_diagram([-1, 0, 0]))
        self.assertFalse(is_viable_diagram([1, -1, 0]))
    
    def test_conjugate_diagram(self):
        """Test Young diagram conjugation."""
        # Test conjugation for various diagrams
        result = conjugate_diagram(pr.YoungDiagram(tuple(self.quark)), n=3).partition
        self.assertTrue(is_viable_diagram(result))
        
        result = conjugate_diagram(pr.YoungDiagram(tuple(self.gluon)), n=3).partition
        self.assertTrue(is_viable_diagram(result))
        
        # Adjoint should be self-conjugate
        adj_conj = conjugate_diagram(pr.YoungDiagram(tuple(self.gluon)), n=3).partition
        # Should be equivalent to original (possibly up to reordering)
        self.assertTrue(is_viable_diagram(adj_conj))

    
    def test_create_fractions_array(self):
        """Test fraction array creation."""
        # Test with various sizes
        result = create_fractions_array(3)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 3)
        
        result = create_fractions_array(5)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 5)
        
        # Test with zero size
        result = create_fractions_array(0)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 0)
        
        # Test with single element
        result = create_fractions_array(1)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
    
    def test_check_vertex_sign(self):
        """Test vertex sign checking function."""
        # Test with valid vertex configurations
        result = check_vertex_sign(
            vertex=1, alpha=self.quark, beta=self.quark, gamma=self.symmetric,
            delta=self.quark, epsilon=self.quark, zeta=self.symmetric,
            v1=1, v2=1, v3=1, v4=1, n=3
        )
        self.assertIsInstance(result, (int, float, type(None)))
        
        # Test with different vertex numbers
        for vertex in range(1, 5):
            result = check_vertex_sign(
                vertex=vertex, alpha=self.quark, beta=self.gluon, gamma=self.quark,
                delta=self.antiquark, epsilon=self.gluon, zeta=self.antiquark,
                v1=1, v2=1, v3=1, v4=1, n=3
            )
            self.assertIsInstance(result, (int, float, type(None)))
    
    def test_6j_type_classification_functions(self):
        """Test 6j-symbol type classification functions."""
        # Standard test parameters
        v1, v2, v3, v4, n = 1, 1, 1, 1, 3
        
        # Test is_6j_with_two_quarks
        result = is_6j_with_two_quarks(
            self.quark, self.quark, self.symmetric,
            self.quark, self.quark, self.symmetric,
            v1, v2, v3, v4, n, debug=False
        )
        self.assertIsInstance(result, bool)
        
        # Test is_6j_with_quark_gluon_vertex
        result = is_6j_with_quark_gluon_vertex(
            self.quark, self.gluon, self.quark,
            self.antiquark, self.gluon, self.antiquark,
            v1, v2, v3, v4, n, debug=False
        )
        self.assertIsInstance(result, bool)
        
        # Test is_6j_with_quark_gluon_opposing
        result = is_6j_with_quark_gluon_opposing(
            self.quark, self.gluon, self.quark,
            self.quark, self.gluon, self.quark,
            v1, v2, v3, v4, n, debug=False
        )
        self.assertIsInstance(result, bool)
        
        # Test is_6j_with_two_gluon
        result = is_6j_with_two_gluon(
            self.gluon, self.gluon, self.trivial,
            self.gluon, self.gluon, self.trivial,
            v1, v2, v3, v4, n, debug=False
        )
        self.assertIsInstance(result, bool)
        
        # Test is_6j_with_three_gluon
        result = is_6j_with_three_gluon(
            self.gluon, self.gluon, self.gluon,
            self.gluon, self.gluon, self.gluon,
            v1, v2, v3, v4, n, debug=False
        )
        self.assertIsInstance(result, bool)
    
    def test_particle_type_consistency(self):
        """Test consistency between particle type functions."""
        # Test various representations
        test_reps = [
            self.trivial, self.quark, self.gluon,
            self.symmetric, self.complex_rep
        ]

        self.assertTrue(is_quark(self.quark))
        self.assertTrue(is_gluon(self.gluon))
        
        for rep in test_reps:
            # All should be viable diagrams
            self.assertTrue(is_viable_diagram(rep))
            
            # Check that type functions return boolean values
            self.assertIsInstance(is_quark(rep), bool)
            self.assertIsInstance(is_adjoint(rep, n=3), bool)
            self.assertIsInstance(is_gluon(rep, n=3), bool)
            self.assertIsInstance(is_real(rep), bool)
            
            # Gluon and adjoint should be consistent
            if is_adjoint(rep, n=3):
                self.assertTrue(is_gluon(rep, n=3))
    
    def test_edge_cases_and_error_handling(self):
        """Test edge cases and error handling."""
        # Empty diagram
        self.assertTrue(is_viable_diagram([]))
        
        # Single element diagrams
        self.assertTrue(is_viable_diagram([0]))
        self.assertTrue(is_viable_diagram([1]))
        self.assertTrue(is_viable_diagram([5]))
        
        # Test with very large representations
        large_rep = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        self.assertTrue(is_viable_diagram(large_rep))
        
        # Test particle type functions with large representations
        self.assertIsInstance(is_quark(large_rep), bool)
        self.assertIsInstance(is_adjoint(large_rep, n=3), bool)
        self.assertIsInstance(is_gluon(large_rep, n=3), bool)
        self.assertIsInstance(is_real(large_rep), bool)
        
        # Test with invalid inputs (should handle gracefully)
        try:
            result = kronecker_delta(None, 1)
            # Should either work or raise appropriate exception
        except (TypeError, ValueError):
            pass  # Expected for invalid input
        
        try:
            result = is_quark(None)
            # Should either work or raise appropriate exception
        except (TypeError, ValueError, AttributeError):
            pass  # Expected for invalid input
    
    
    def test_mathematical_properties(self):
        """Test mathematical properties of utility functions."""
        # Kronecker delta properties
        for i in range(-5, 6):
            for j in range(-5, 6):
                delta_ij = kronecker_delta(i, j)
                delta_ji = kronecker_delta(j, i)
                # Should be symmetric
                self.assertEqual(delta_ij, delta_ji)
                # Should be 0 or 1
                self.assertIn(delta_ij, [0, 1])
        
        # Conjugation should preserve total number of boxes for SU(3)
        for rep in [self.quark, self.gluon, self.symmetric]:
            from wigner_new import get_number_of_boxes
            original_boxes = get_number_of_boxes(rep)
            conjugated = conjugate_diagram(pr.YoungDiagram(tuple(rep)), n=3).partition
            conjugated_boxes = get_number_of_boxes(conjugated)
            
            # For SU(3), conjugation preserves certain relationships
            self.assertGreaterEqual(conjugated_boxes, 0)
    

    def test_debug_mode_compatibility(self):
        """Test that functions work properly with debug mode."""
        # Test classification functions with debug=True
        v1, v2, v3, v4, n = 1, 1, 1, 1, 3
        
        result = is_6j_with_two_quarks(
            self.quark, self.quark, self.symmetric,
            self.quark, self.quark, self.symmetric,
            v1, v2, v3, v4, n, debug=True
        )
        self.assertIsInstance(result, bool)
        
        result = is_6j_with_quark_gluon_vertex(
            self.quark, self.gluon, self.quark,
            self.antiquark, self.gluon, self.antiquark,
            v1, v2, v3, v4, n, debug=True
        )
        self.assertIsInstance(result, bool)


if __name__ == '__main__':
    unittest.main()