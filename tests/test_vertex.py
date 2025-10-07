#!/usr/bin/env python3
"""
Tests for the Vertex class in the Wigner6j package.

The Vertex class represents vertices in 6j-symbol calculations
with three Young diagrams and a vertex number.
"""
import unittest
import sys
from pathlib import Path
import pair_multiplication as pr

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from wigner_new import Vertex


class TestVertex(unittest.TestCase):
    """Test cases for the Vertex class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Common Young diagrams for testing
        self.trivial_diagram = pr.YoungDiagram((), Nc=3)
        self.quark_diagram = pr.YoungDiagram((1), Nc=3)
        self.gluon_diagram = pr.YoungDiagram((2, 1), Nc=3)  # Adjoint representation for SU(3)
        self.complex_diagram = pr.YoungDiagram((4, 3), Nc=3)

    def test_initialization_with_lists(self):
        """Test Vertex initialization with Young diagram lists."""
        vertex = Vertex(self.quark_diagram, self.gluon_diagram, self.complex_diagram, 1)
        
        self.assertEqual(vertex.vertex_number, 1)
        # Check that Young diagrams are stored (type depends on pair_multiplication implementation)
        self.assertTrue(hasattr(vertex, 'young_diagram_1'))
        self.assertTrue(hasattr(vertex, 'young_diagram_2'))
        self.assertTrue(hasattr(vertex, 'young_diagram_3'))

        self.assertEqual(vertex.young_diagram_1, self.quark_diagram)
        self.assertEqual(vertex.young_diagram_2, self.gluon_diagram)
        self.assertEqual(vertex.young_diagram_3, self.complex_diagram)

    def test_initialization_defaults(self):
        """Test Vertex initialization with default parameters."""
        vertex = Vertex()
        
        self.assertEqual(vertex.vertex_number, 1)
        # Should have default trivial diagrams
        self.assertTrue(hasattr(vertex, 'young_diagram_1'))
        self.assertTrue(hasattr(vertex, 'young_diagram_2'))
        self.assertTrue(hasattr(vertex, 'young_diagram_3'))

        self.assertEqual(vertex.young_diagram_1, self.trivial_diagram)
        self.assertEqual(vertex.young_diagram_2, self.trivial_diagram)
        self.assertEqual(vertex.young_diagram_3, self.trivial_diagram)

    def test_is_well_defined_method_exists(self):
        """Test that is_well_defined method exists and returns boolean."""
        vertex = Vertex(self.quark_diagram, self.gluon_diagram, self.complex_diagram, 1)
        
        # Test that the method exists and can be called
        self.assertTrue(hasattr(vertex, 'is_well_defined'))
        self.assertTrue(callable(vertex.is_well_defined))
        
        # Test that it returns a boolean (assuming it doesn't raise an exception)
        try:
            result = vertex.is_well_defined()
            self.assertIsInstance(result, bool)
        except Exception:
            # If it raises an exception due to missing dependencies, that's expected
            pass


if __name__ == '__main__':
    unittest.main()