r"""
Unit tests for the CLI.
"""

import os
import tempfile
import unittest


class TestCLI(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory(dir='/tmp')
        self.cif_dirname = 'tests/CIFs'

    def test_cli(self):
        os.system(f'moxel create -g 5 {self.cif_dirname} {self.tempdir.name}')
        os.system(f'moxel clean {self.tempdir.name}')
        os.system(f'moxel prepare {self.tempdir.name}/clean_names.json')

        # Check that the files are correctly created.
        os.path.isfile(f'{self.tempdir.name}/clean_voxels.npy')
        os.path.isfile(f'{self.tempdir.name}/clean_names.json')
        for mode in ['train', 'validation', 'test']:
            self.assertTrue(os.path.isfile(f'{self.tempdir.name}/{mode}.json'))

        # Check that LightningCLI works.
        ...

    def tearDown(self):
        self.tempdir.cleanup()


if __name__ == '__main__':
    unittest.main()
