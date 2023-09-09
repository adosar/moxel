# This file is part of MOXελ.
# Copyright (C) 2023 Antonios P. Sarikas

# MOXελ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
Run tests from the project's root.

python -m unittest tests.<test_module>
"""

import os
import sys
import json
import unittest
import tempfile
import numpy as np
from pymatgen.core import Structure
sys.path.insert(0, './src')
from moxel import *


class TestMoxelUtils(unittest.TestCase):
    
    def test_Grid(self):
        cif_pathname = 'tests/CIFs/foo/IRMOF-1.cif'
        grid_size = 5
        epsilon = 40
        sigma = 5
        cutoff = 12

        grid = Grid(grid_size=grid_size, epsilon=epsilon, cutoff=cutoff, sigma=sigma)
        grid.load_structure(cif_pathname)
        grid.calculate(cubic_box=True)

        # Check that the name of the structure is correct.
        self.assertEqual(grid.structure_name, 'IRMOF-1')

        # Check that the voxels have the correct shape.
        self.assertEqual(grid.voxels.shape, (grid_size,)*3)

        # Check that attributes are correctly set.
        self.assertEqual(grid.grid_size, grid_size)
        self.assertEqual(grid.epsilon, epsilon)
        self.assertEqual(grid.sigma, sigma)
        self.assertEqual(grid.cutoff, cutoff)
        self.assertTrue(grid.cubic_box)

    def test_mic_scale_factors(self):
        # Load a structure that must be scaled.
        cif_pathname = 'tests/CIFs/foo/ZnMOF-74.cif'
        structure = Structure.from_file(cif_pathname)
        scale_factors = mic_scale_factors(10, structure.lattice.matrix)

        # At least one scale factor should be different than 1.
        self.assertTrue(np.any(scale_factors != 1))

        # Load a structure that must not be scaled.
        cif_pathname = 'tests/CIFs/foo/IRMOF-1.cif'
        structure = Structure.from_file(cif_pathname)
        scale_factors = mic_scale_factors(10, structure.lattice.matrix)

        # At least one scale factor should be different than 1.
        self.assertTrue(np.all(scale_factors == 1))

    def test_voxels_from_file(self):
        cif_pathname = 'tests/CIFs/foo/IRMOF-1.cif'
        grid_size = 10
        cutoff = 9
        epsilon = 49
        sigma = 2

        # `out_1` -> `cubic_box=True` & `length=15`.
        out_1 = voxels_from_file(
                cif_pathname, grid_size=grid_size,
                cutoff=cutoff, epsilon=epsilon,
                sigma=sigma, cubic_box=True,
                length=30, only_voxels=False
                )

        # Check that arguments are properly passed.
        self.assertTrue(isinstance(out_1, Grid))
        self.assertEqual(out_1.grid_size, grid_size)
        self.assertEqual(out_1.cutoff, cutoff)
        self.assertEqual(out_1.epsilon, epsilon)
        self.assertEqual(out_1.sigma, sigma)
        self.assertTrue(out_1.cubic_box)
        self.assertTrue(np.all(out_1.voxels.shape ==  np.array([grid_size]*3)))

        # `out_2` -> `cubic_box=False`.
        out_2 = voxels_from_file(
                cif_pathname, grid_size=grid_size,
                cutoff=cutoff, epsilon=epsilon,
                sigma=sigma, cubic_box=False,
                only_voxels=False
                )
        # Check that `cubic_box` works.
        self.assertFalse(np.all(out_1.voxels == out_2.voxels))

        # `out_3` -> `cubix_box=True` & `length=10`.
        out_3 = voxels_from_file(
                cif_pathname, grid_size=grid_size,
                cutoff=cutoff, epsilon=epsilon,
                sigma=sigma, cubic_box=True,
                length=10, only_voxels=False
                )
        # Check that `length` works.
        self.assertFalse(np.all(out_1.voxels == out_3.voxels))


    def test_voxels_from_files(self):
        n_files = 3
        grid_size = 10
        cif_dir = 'tests/CIFs/foo'
        cif_files = [f'{cif_dir}/{i}' for i in os.listdir(cif_dir)][:n_files]

        # Check that output file is properly stored.
        with tempfile.TemporaryDirectory() as dir_path:
            out_pathname = f'{dir_path}/foo.npy'
            voxels_from_files(cif_files, out_pathname=out_pathname, grid_size=grid_size)

            voxels = np.load(out_pathname, mmap_mode='r')

            self.assertTrue(
                    np.all(voxels.shape == np.array([n_files, *[grid_size]*3]))
                    )

    def test_voxels_from_dir(self):
        # The `cif_dir` contains 1 corrupted .cif file.
        cif_dir = 'tests/CIFs/foo'
        grid_size = 10

        with tempfile.TemporaryDirectory() as dir_path:
            out_pathname = f'{dir_path}/foo.npy'
            voxels_from_dir(cif_dir, grid_size=grid_size, out_pathname=out_pathname)

            voxels = np.load(out_pathname, mmap_mode='r')

            # Check that corrupted .cif files aren't processed.
            bad_cifs = [np.all(x == 0) for x in voxels]
            self.assertEqual(
                    sum(bad_cifs),
                    len([i for i in os.listdir(cif_dir) if i.startswith('corrupted')])
                    )

    def test_batch_clean_and_merge_single(self):
        cif_dir = 'tests/CIFs/bar'
        grid_size = 10
        names = sorted(os.listdir(cif_dir))
        clean_names_init = [i for i in names if not i.startswith('corrupted')]

        with tempfile.TemporaryDirectory() as dir_path:
            out_pathname = f'{dir_path}/voxels.npy'
            voxels_from_dir(cif_dir, grid_size=grid_size, out_pathname=out_pathname)
            voxels_shape = np.array([len(names), *[grid_size]*3])

            # Add the names under the batch directory.
            with open(f'{dir_path}/names.json', 'w') as fhand:
                json.dump({'names': names}, fhand)

            exit_status = batch_clean_and_merge([dir_path])

            # Check the exit status.
            self.assertEqual(exit_status, 1)

            clean_voxels = np.load(f'{dir_path}/clean_voxels.npy', mmap_mode='r')

            # Load the clean names.
            with open(f'{dir_path}/clean_names.json', 'r') as fhand:
                clean_names_final = json.load(fhand)['names']

            # Check that corrupted voxels/names are deleted .
            self.assertFalse(voxels_shape[0] == clean_voxels.shape[0])
            self.assertTrue(np.all(voxels_shape[1:] == clean_voxels.shape[1:]))
            self.assertEqual(clean_names_init, clean_names_final)

            # Check that the remaining voxels are properly processed.
            self.assertEqual(sum([i for i in clean_voxels if np.all(i == 0)]), 0)

            # Check that no further processing is required.
            os.remove(f'{dir_path}/names.json')
            os.remove(f'{dir_path}/voxels.npy')
            os.rename(f'{dir_path}/clean_voxels.npy', f'{dir_path}/voxels.npy')
            os.rename(f'{dir_path}/clean_names.json', f'{dir_path}/names.json')

            exit_status = batch_clean_and_merge([dir_path])

            self.assertEqual(exit_status, 0)

    def test_batch_clean_and_merge_multiple(self):
        grid_size = 10
        cif_dir_1 = 'tests/CIFs/foo'
        cif_dir_2 = 'tests/CIFs/bar'

        names_1 = sorted(os.listdir(cif_dir_1))
        names_2 = sorted(os.listdir(cif_dir_2))

        names = names_1 + names_2
        names_corrupted = [i for i in names if i.startswith('corrupted')]

        # Order and repetition matters (avoid using `set`)
        names_clean_init = [j for j in names if j not in names_corrupted]

        with tempfile.TemporaryDirectory() as dir_path:

            # Creating the 1st batch.
            tempdir_1 = tempfile.TemporaryDirectory()
            dir_path_1 = tempdir_1.name
            out_pathname_1 = f'{dir_path_1}/voxels.npy'
            voxels_from_dir(cif_dir_1, grid_size=grid_size, out_pathname=out_pathname_1)

            with open(f'{dir_path_1}/names.json', 'w') as fhand:
                json.dump({'names': names_1}, fhand)

            # Creating the 2nd batch.
            tempdir_2 = tempfile.TemporaryDirectory()
            dir_path_2 = tempdir_2.name
            out_pathname_2 = f'{dir_path_2}/voxels.npy'
            voxels_from_dir(cif_dir_2, grid_size=grid_size, out_pathname=out_pathname_2)

            with open(f'{dir_path_2}/names.json', 'w') as fhand:
                json.dump({'names': names_2}, fhand)

            exit_status = batch_clean_and_merge([dir_path_1, dir_path_2], out_pathname=dir_path)

            # Check that processing is required.
            self.assertEqual(exit_status, 1)

            # Check that no further processing is required.
            os.rename(f'{dir_path}/clean_voxels.npy', f'{dir_path}/voxels.npy')
            os.rename(f'{dir_path}/clean_names.json', f'{dir_path}/names.json')

            exit_status = batch_clean_and_merge([dir_path])
            self.assertEqual(exit_status, 0)

            # Check that the correct names are stored.
            with open(f'{dir_path}/clean_names.json', 'r') as fhand:
                names_clean_final = json.load(fhand)['names']

            self.assertEqual(names_clean_init, names_clean_final)


if __name__ == '__main__':
    unittest.main()
