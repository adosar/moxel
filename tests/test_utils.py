import os
import sys
import json
import unittest
import tempfile
import numpy as np
from pymatgen.core import Structure
sys.path.insert(0, '/home/asar/projects/moxel_pypi/src')
from moxel import *


class TestMoxelUtils(unittest.TestCase):
    
    def test_mic_scale_factors(self):
        # Load a structure that must be scaled
        cif_pathname = 'tests/CIFs/foo/ZnMOF-74.cif'
        structure = Structure.from_file(cif_pathname)
        scale_factors = mic_scale_factors(10, structure.lattice.matrix)

        # At least one scale factor should be different than 1
        self.assertTrue(np.any(scale_factors != 1))

        # Load a structure that must not be scaled
        cif_pathname = 'tests/CIFs/foo/IRMOF-1.cif'
        structure = Structure.from_file(cif_pathname)
        scale_factors = mic_scale_factors(10, structure.lattice.matrix)

        # At least one scale factor should be different than 1
        self.assertTrue(np.all(scale_factors == 1))

    def test_voxels_from_file(self):
        cif_pathname = 'tests/CIFs/foo/IRMOF-1.cif'
        grid_size = 10
        cutoff = 9
        epsilon = 49
        sigma = 2

        # `out_1` -> `cubic_box=True` & `length=15`
        out_1 = voxels_from_file(
                cif_pathname, grid_size=grid_size,
                cutoff=cutoff, epsilon=epsilon,
                sigma=sigma, cubic_box=True,
                length=30, only_voxels=False
                )
        # Check that arguments are properly passed
        self.assertTrue(isinstance(out_1, Grid))
        self.assertEqual(out_1.grid_size, grid_size)
        self.assertEqual(out_1.cutoff, cutoff)
        self.assertEqual(out_1.epsilon, epsilon)
        self.assertEqual(out_1.sigma, sigma)
        self.assertTrue(out_1.cubic_box)
        self.assertTrue(np.all(out_1.voxels.shape ==  np.array([grid_size]*3)))

        # `out_2` -> `cubic_box=False`
        out_2 = voxels_from_file(
                cif_pathname, grid_size=grid_size,
                cutoff=cutoff, epsilon=epsilon,
                sigma=sigma, cubic_box=False,
                only_voxels=False
                )
        # Check that `cubic_box` works
        self.assertFalse(np.all(out_1.voxels == out_2.voxels))

        # `out_3` -> `cubix_box=True` & `length=10`
        out_3 = voxels_from_file(
                cif_pathname, grid_size=grid_size,
                cutoff=cutoff, epsilon=epsilon,
                sigma=sigma, cubic_box=True,
                length=10, only_voxels=False
                )
        # Check that `length` works
        self.assertFalse(np.all(out_1.voxels == out_3.voxels))


    def test_voxels_from_files(self):
        n_files = 3
        grid_size = 10
        cif_dir = 'tests/CIFs/foo'
        cif_files = [f'{cif_dir}/{i}' for i in os.listdir(cif_dir)][:n_files]

        # Check that output file is properly stored
        with tempfile.TemporaryDirectory() as dir_path:
            out_name = f'{dir_path}/foo.npy'
            voxels_from_files(cif_files, out_name=out_name, grid_size=grid_size)

            voxels = np.load(out_name, mmap_mode='r')

            self.assertTrue(
                    np.all(voxels.shape == np.array([n_files, *[grid_size]*3]))
                    )

    def test_voxels_from_dir(self):
        # The `cif_dir` contains 1 corrupted .cif file
        cif_dir = 'tests/CIFs/foo'
        grid_size = 10

        with tempfile.TemporaryDirectory() as dir_path:
            out_name = f'{dir_path}/foo.npy'
            voxels_from_dir(cif_dir, grid_size=grid_size, out_name=out_name)

            voxels = np.load(out_name, mmap_mode='r')

            # Check that corrupted .cif files aren't processed
            bad_cifs = [np.all(x == 0) for x in voxels]
            self.assertEqual(
                    sum(bad_cifs),
                    len([i for i in os.listdir(cif_dir) if i.startswith('corrupted')])
                    )

    def test_batch_clean_and_merge_single(self):
        cif_dir = 'tests/CIFs/foo'
        grid_size = 10
        cif_names = os.listdir(cif_dir)
        n_corrupted = len([i for i in cif_names if i.startswith('corrupted')])


        with tempfile.TemporaryDirectory() as dir_path:
            out_name = f'{dir_path}/voxels.npy'
            voxels_from_dir(cif_dir, grid_size=grid_size, out_name=out_name)
            voxels_shape = np.array([len(cif_names), *[grid_size]*3])

            # Add the names under the batch directory
            with open(f'{dir_path}/names.json', 'w') as fhand:
                json.dump({'names': cif_names}, fhand)

            exit_status = batch_clean_and_merge([dir_path])

            # Check the exit status
            self.assertEqual(exit_status, 1)

            clean_voxels = np.load(f'{dir_path}/clean_voxels.npy', mmap_mode='r')

            # Load the clean names
            with open(f'{dir_path}/clean_names.json', 'r') as fhand:
                clean_names = json.load(fhand)['names']

            # Check that corrupted voxels/names are deleted 
            self.assertFalse(voxels_shape[0] == clean_voxels.shape[0])
            self.assertTrue(np.all(voxels_shape[1:] == clean_voxels.shape[1:]))
            self.assertEqual(len(cif_names) - n_corrupted, len(clean_names))

            # Check that the remaining voxels are properly processed
            self.assertEqual(sum([i for i in clean_voxels if np.all(i == 0)]), 0)

    def test_batch_clean_and_merge_multiple
        cif_dir_1 = 'tests/CIFs/foo'
        cif_dir_2 = 'tests/CIFs/bar'
        cif_names = os.listdir(cif_dir_1) + os.listdir(cif_dir_2)
        n_corrupted = [len([i for i in cif_names if i.startswith('corrupted')]

        with tempfile.TemporaryDirectory() as dir_path:
            # Creating the 1st batch
            with tempfile.TemporaryDirectory() as dir_path_1
                voxels_from_dir(cif_dir_1


if __name__ == '__main__':
    unittest.main()

