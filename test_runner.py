#!/usr/bin/env python3
"""
Test runner for Wigner6j package tests.

Run all tests:
    python test_runner.py

Run specific test modules:
    python test_runner.py test_fraction
    python test_runner.py test_vertex test_wigner
"""
import sys
import unittest
import os
from pathlib import Path

# Add the parent directory to Python path so we can import the modules
sys.path.insert(0, str(Path(__file__).parent.parent))

def run_tests(test_modules=None):
    """Run the test suite."""
    # Discover tests in the tests directory
    loader = unittest.TestLoader()
    
    if test_modules:
        # Run specific test modules
        suite = unittest.TestSuite()
        for module_name in test_modules:
            try:
                module = __import__(f'tests.{module_name}', fromlist=[''])
                suite.addTests(loader.loadTestsFromModule(module))
            except ImportError as e:
                print(f"Warning: Could not import test module '{module_name}': {e}")
    else:
        # Run all tests
        suite = loader.discover('tests', pattern='test_*.py')
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()

if __name__ == '__main__':
    test_modules = sys.argv[1:] if len(sys.argv) > 1 else None
    success = run_tests(test_modules)
    sys.exit(0 if success else 1)