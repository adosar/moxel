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
import unittest
import tempfile
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from moxel.utils import (
        Grid, mic_scale_factors, voxels_from_file,
        voxels_from_files, voxels_from_dir
        )


class TestUtils(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory(dir='/tmp')

    def test_grid(self):
        cif_pathname = 'tests/CIFs/IRMOF-1.cif'
        grid_size = 5
        epsilon = 40
        sigma = 1e4  # For repulsive energies to go inf.
        cutoff = 12
        cubic_box = 35

        grid = Grid(grid_size=grid_size, epsilon=epsilon, cutoff=cutoff, sigma=sigma)
        grid.load_structure(cif_pathname)
        grid.calculate(cubic_box=cubic_box)

        # Check that the name of the structure is correct.
        self.assertEqual(grid.structure_name, 'IRMOF-1')

        # Check that the voxels have the correct shape.
        self.assertEqual(grid.voxels.shape, (grid_size,)*3)

        # Check that attributes are correctly set.
        self.assertEqual(grid.grid_size, grid_size)
        self.assertEqual(grid.epsilon, epsilon)
        self.assertEqual(grid.sigma, sigma)
        self.assertEqual(grid.cutoff, cutoff)
        self.assertEqual(grid.cubic_box, cubic_box)

        # Check that values are filled properly.
        self.assertEqual(grid.voxels.max(), np.inf)

    def test_mic_scale_factors(self):
        # Load a structure that must be scaled.
        cif_pathname = 'tests/CIFs/ZnMOF-74.cif'
        structure = Structure.from_file(cif_pathname)
        scale_factors = mic_scale_factors(10, structure.lattice.matrix)

        # At least one scale factor should be different than 1.
        self.assertTrue(np.any(scale_factors != 1))

        # Load a structure that shouldn't be scaled.
        cif_pathname = 'tests/CIFs/IRMOF-1.cif'
        structure = Structure.from_file(cif_pathname)
        scale_factors = mic_scale_factors(10, structure.lattice.matrix)

        # All scale factors should be 1.
        self.assertTrue(np.all(scale_factors == 1))

    def test_voxels_from_file(self):
        cif_pathname = 'tests/CIFs/IRMOF-1.cif'
        grid_size = 5
        cutoff = 9
        epsilon = 49
        sigma = 2
        cubic_box = 25

        out = voxels_from_file(
                cif_pathname, grid_size=grid_size,
                cutoff=cutoff, epsilon=epsilon,
                sigma=sigma, cubic_box=cubic_box,
                only_voxels=False
                )

        # Check that arguments are properly passed.
        self.assertTrue(isinstance(out, Grid))
        self.assertEqual(out.grid_size, grid_size)
        self.assertEqual(out.cutoff, cutoff)
        self.assertEqual(out.epsilon, epsilon)
        self.assertEqual(out.sigma, sigma)
        self.assertEqual(out.cubic_box, cubic_box)

        # Check that output shape is correct.
        self.assertEqual(out.voxels.shape, (grid_size,)*3)

    def test_voxels_from_files(self):
        # The test assumes all files are processable.
        grid_size = 5
        cif_dir = 'tests/CIFs'
        cif_files = [f'{cif_dir}/{i}' for i in os.listdir(cif_dir)]
        names = [Path(i).stem for i in cif_files]

        # Calculate voxels and store them.
        out_pathname = f'{self.tempdir.name}/voxels_data'
        voxels_from_files(cif_files, out_pathname, grid_size=grid_size)

        npy_files = [os.path.join(out_pathname, f) for f in os.listdir(out_pathname)]

        # Check that all files are processed properly.
        for f in npy_files:
            voxels  = np.load(f)
            self.assertTrue(Path(f).stem in names)
            self.assertTrue(voxels.shape, (grid_size,)*3)

    def test_voxels_from_dir(self):
        # The test assumes all files are processable.
        cif_dir = 'tests/CIFs'
        grid_size = 5

        # Calculate voxels and store them.
        out_pathname = f'{self.tempdir.name}/voxels_data'
        voxels_from_dir(cif_dir, out_pathname, grid_size=grid_size)

        # Check that all files are processed.
        self.assertEqual(len(os.listdir(out_pathname)), len(os.listdir(cif_dir)))

    def test_load_file(self):
        cif_dir = 'tests/CIFs'

        grid = Grid()
        grid.load_structure(f'{cif_dir}/IRMOF-1.cif')

        grid = Grid()
        grid.load_structure(f'{cif_dir}/MnH28C26(N2Cl)2.json')

    def tearDown(self):
        self.tempdir.cleanup()


if __name__ == '__main__':
    unittest.main()
