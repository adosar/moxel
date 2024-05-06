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

import os
import sys
import json
import unittest
import tempfile
from pathlib import Path
import numpy as np
from pymatgen.core import Structure
from moxel.utils import *


class TestUtils(unittest.TestCase):
    def test_grid(self):
        cif_pathname = 'tests/CIFs/IRMOF-1.cif'
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

    def test_mic_scale_factors(self):
        # Load a structure that must be scaled.
        cif_pathname = 'tests/CIFs/ZnMOF-74.cif'
        structure = Structure.from_file(cif_pathname)
        scale_factors = mic_scale_factors(10, structure.lattice.matrix)

        # At least one scale factor should be different than 1.
        self.assertTrue(np.any(scale_factors != 1))

        # Load a structure that must not be scaled.
        cif_pathname = 'tests/CIFs/IRMOF-1.cif'
        structure = Structure.from_file(cif_pathname)
        scale_factors = mic_scale_factors(10, structure.lattice.matrix)

        # At least one scale factor should be different than 1.
        self.assertTrue(np.all(scale_factors == 1))

    def test_voxels_from_file(self):
        cif_pathname = 'tests/CIFs/IRMOF-1.cif'
        grid_size = 5
        cutoff = 9
        epsilon = 49
        sigma = 2
        cubic_box = True
        length = 22

        out = voxels_from_file(
                cif_pathname, grid_size=grid_size,
                cutoff=cutoff, epsilon=epsilon,
                sigma=sigma, cubic_box=cubic_box,
                length=length, only_voxels=False
                )

        # Check that arguments are properly passed.
        self.assertTrue(isinstance(out, Grid))
        self.assertEqual(out.grid_size, grid_size)
        self.assertEqual(out.cutoff, cutoff)
        self.assertEqual(out.epsilon, epsilon)
        self.assertEqual(out.sigma, sigma)

        # Check that output is correct.
        self.assertEqual(out.voxels.shape, (grid_size,)*3)
        self.assertFalse(np.all(out.voxels == 0))

    def test_voxels_from_files(self):
        grid_size = 5
        cif_dir = 'tests/CIFs'
        cif_files = [f'{cif_dir}/{i}' for i in os.listdir(cif_dir)]
        n_files = len(cif_files)

        names = [Path(i).stem for i in cif_files]

        # Check that output files are properly stored.
        with tempfile.TemporaryDirectory() as dir_path:
            voxels_from_files(cif_files, out_pathname=dir_path, grid_size=grid_size)

            voxels = np.load(f'{dir_path}/voxels.npy', mmap_mode='r')
            stored_names = get_names(f'{dir_path}/names.json')

            self.assertEqual(voxels.shape, (n_files, grid_size, grid_size, grid_size))
            self.assertEqual(names, stored_names)

    def test_voxels_from_dir(self):
        # The cif_dir contains 2 corrupted .cif files.
        cif_dir = 'tests/CIFs'
        grid_size = 5

        cif_names = sorted(os.listdir(cif_dir))
        n_files = len(cif_names)
        n_bad_files = len([i for i in cif_names if i.startswith('corrupted')])

        names = [Path(i).stem for i in cif_names]

        with tempfile.TemporaryDirectory() as dir_path:
            voxels_from_dir(cif_dir, dir_path, grid_size=grid_size)

            voxels = np.load(f'{dir_path}/voxels.npy', mmap_mode='r')
            stored_names = get_names(f'{dir_path}/names.json')

            self.assertEqual(voxels.shape, (n_files, grid_size, grid_size, grid_size))
            self.assertEqual(names, stored_names)

            # Check that corrupted .cif files aren't processed.
            n_bad_voxels = sum([np.all(x == 0) for x in voxels])
            self.assertEqual(n_bad_files, n_bad_voxels)

    def test_batch_clean(self):
        cif_dir = 'tests/CIFs'
        grid_size = 5

        with tempfile.TemporaryDirectory() as dir_path:
            voxels_from_dir(cif_dir, dir_path, grid_size=grid_size)
            
            batch_clean(dir_path)

            voxels = np.load(f'{dir_path}/voxels.npy', mmap_mode='r')
            names = get_names(f'{dir_path}/names.json')

            clean_voxels = np.load(f'{dir_path}/clean_voxels.npy', mmap_mode='r')
            clean_names = get_names(f'{dir_path}/clean_names.json')

            # Get the indices of filled voxels.
            filled_idx = [i for i, x in enumerate(voxels) if not np.all(x == 0)]

            # Check that the voxels have been cleaned.
            self.assertTrue(np.all(voxels[filled_idx] == clean_voxels))

            # Check that the names have been correctly stored.
            self.assertEqual([names[i] for i in filled_idx], clean_names)

    def test_load_file(self):
        cif_dir = 'tests/CIFs'

        grid = Grid()
        grid.load_structure(f'{cif_dir}/IRMOF-1.cif')

        grid = Grid()
        grid.load_structure(f'{cif_dir}/MnH28C26(N2Cl)2.json')


if __name__ == '__main__':
    unittest.main()
