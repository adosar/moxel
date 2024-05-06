import unittest
import doctest
import moxel.visualize


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(moxel.visualize))
    return tests


if __name__ == '__main__':
    unittest.main()
