"""
Run tests from the project's root.

python -m unittest tests.<test_module>
"""

import os
import unittest
import tempfile


class TestMoxelCLI(unittest.TestCase):

    def test_help(self):
        exit_status = os.system('python -m src.moxel -h')
        self.assertEqual(exit_status, 0)

    def test_run(self):
        with tempfile.TemporaryDirectory() as dir_path:
            file_name = 'voxels.npy'
            exit_status = os.system(
                    f'python -m src.moxel tests/CIFs/foo\
                            -n 5 -o {dir_path}/{file_name}'
                    )
            self.assertEqual(exit_status, 0)


if __name__ == '__main__':
    unittest.main()
