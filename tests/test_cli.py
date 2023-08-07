import os
import unittest
import tempfile


class TestMoxelHelp(unittest.TestCase):

    def test_help(self):
        exit_status = os.system('python -m src.moxel -h')
        self.assertEqual(exit_status, 0)


class TestMoxelRun(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def test_run(self):
        exit_status = os.system(
                f'python -m src.moxel tests/cif_directory\
                        -o {self.tempdir}/voxels.npy'
                )
        self.assertEqual(exit_status, 0)
        #self.assertEqual(os.system(f'ls 

    #def tearDown(self):
    #    os
    #    pass


if __name__ == '__main__':
    unittest.main()
