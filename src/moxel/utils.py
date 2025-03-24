# This file is part of MOXελ.
# Copyright (C) 2023-2024 Antonios P. Sarikas

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

r"""
This module provides helper functions for creating voxels.

.. note::
    Currently, interactions are modelled with the Lennard-Jones (LJ) potential.

.. attention::
    Consider tuning the ``n_jobs`` parameter to get the best performance
    for your system::

        from timeit import timeit

        setup = 'from moxel.utils import voxels_from_file'
        n_jobs = [1, 2, 8, 16]  # Modify this according to your system.

        for n in n_jobs:
            stmt = f'voxels_from_file("path/to/cif", n_jobs={n})'
            time = timeit(stmt=stmt, setup=setup, number=1)
            print(f'Time with {n} jobs: {time:.3f} s')
"""

import os
import itertools
from pathlib import Path
from multiprocessing import Pool
import warnings

import numpy as np
from tqdm import tqdm
from pymatgen.core import Structure

from . _params import lj_params

warnings.filterwarnings('ignore')

# Default values for voxels calculation.
GRID_SIZE = 25
CUTOFF = 10.
EPSILON = 50.
SIGMA = 2.5
CUBIC_BOX = None
N_JOBS = None


def mic_scale_factors(r, lattice_vectors):
    r"""
    Return scale factors to satisfy minimum image convention [MIC]_.

    Parameters
    ----------
    r : float
        The cutoff radius used in MIC convetion.
    lattice_vectors : array of shape (3, 3)
        The lattice vectors of the unit cell.
        Each row corresponds to a lattice vector.

    Returns
    -------
    scale_factors : array of shape (3,)
        ``scale_factors[i]`` scales ``lattice_vectors[i]``.

    References
    ----------
    .. [MIC] W. Smith, "The Minimum Image Convention in Non-Cubic MD Cells", 1989.
    """
    a, b, c = lattice_vectors
    volume = np.linalg.norm(np.dot(a, np.cross(b, c)))

    w_a = volume/np.linalg.norm(np.cross(b, c))
    w_b = volume/np.linalg.norm(np.cross(a, c))
    w_c = volume/np.linalg.norm(np.cross(a, b))

    return np.ceil(2 * r / np.array([w_a, w_b, w_c]))


class Grid:
    r"""
    A 3D energy grid over a crystal structure.

    Parameters
    ----------
    grid_size : int, default=25
        Number of grid points along each dimension.
    cutoff : float, default=10.0
        Cutoff radius (Å) for the LJ potential.
    epsilon : float, default=50.0
        Epsilon value (ε/K) of the probe atom.
    sigma : float, default=2.5
        Sigma value (σ/Å) of the probe atom.

    Attributes
    ----------
    structure : :class:`pymatgen.core.structure.Structure`
        Available only after :meth:`Grid.load_structure` has been called.
    structure_name : str
        Available only after :meth:`Grid.load_structure` has been called.
    cubic_box : float or None
        Available only after :meth:`Grid.calculate` has been called.
    voxels : array of shape (grid_size,)*3
       Available only after :meth:`Grid.calculate` has been called.
    """
    def __init__(
            self,
            grid_size=GRID_SIZE,
            *,
            cutoff=CUTOFF,
            epsilon=EPSILON,
            sigma=SIGMA
            ):
        self.grid_size = grid_size
        self.cutoff = cutoff
        self.epsilon = epsilon
        self.sigma = sigma

    def load_structure(self, pathname):
        r"""
        Load a crystal structure from a file in a format supported by
        :meth:`pymatgen.core.Structure.from_file`.

        Parameters
        ----------
        pathname : str
           Pathname to the file.
        """
        self.structure = Structure.from_file(pathname)
        self.structure_name = Path(pathname).stem

    def calculate(
            self,
            cubic_box=CUBIC_BOX,
            potential='lj',
            n_jobs=N_JOBS,
            ):
        r"""
        Iterate over the grid and return voxels.

        For computational efficiency and to assure (approximately) the same
        spatial resolution, the grid is overlayed over a supercell scaled
        according to MIC, see :func:`mic_scale_factors`.

        If lattice angles are significantly different than 90°, to avoid
        distortions set ``cubic_box``. In this case, the grid is overlayed over
        a cubic box of size ``cubic_box`` centered at the origin but periodicity
        is no longer guaranteed.

        Parameters
        ----------
        potential : str, default='lj'
            The potential used to calculate voxels. Currently, only the
            LJ potential is supported.
        cubic_box : float or None, default=None
            If ``None``, the simulation box is a supercell scaled according to
            MIC. Otherwise, cubic box of size ``cubic_box``.
        n_jobs : int, optional
            Number of jobs to run in parallel. If ``None``, then the number returned
            by ``os.cpu_count()`` is used.

        Returns
        -------
        voxels : array of shape (grid_size,)*3
        """
        self.cubic_box = cubic_box

        if cubic_box is not None:
            d = cubic_box / 2
            probe_coords = np.linspace(0 - d, 0 + d, self.grid_size, endpoint=False)  # Cartesian
            self._simulation_box = self.structure
        else:
            probe_coords = np.linspace(0, 1, self.grid_size, endpoint=False)  # Fractional
            scale = mic_scale_factors(self.cutoff, self.structure.lattice.matrix)
            self._simulation_box = self.structure * scale

        if potential == 'lj':
            # Cache LJ parameters for all atoms in the simulation box.
            self._lj_params = np.array(
                    [lj_params[atom.species_string] for atom in self._simulation_box]
                    )

            # Cache fractional coordinates since this is a slow function in pymatgen.
            self._frac_coords = self._simulation_box.frac_coords

            # Embarrassingly parallel.
            with Pool(processes=n_jobs) as p:
                energies = p.map(
                        self.lj_potential, itertools.product(*(probe_coords,)*3)
                        )

        self.voxels = np.array(energies, dtype=np.float32).reshape((self.grid_size,)*3)

        return self.voxels

    def lj_potential(self, coords):
        r"""
        Calculate LJ potential at cartesian or fractional
        coordinates.

        Parameters
        ----------
        coordinates : array_like of shape (3,)
            If ``cubic_box=None`` fractional, else cartesian.

        Returns
        -------
        energy : float
        """
        if self.cubic_box is not None:
            cartesian_coords = coords
        else:
            cartesian_coords = self._simulation_box._lattice.get_cartesian_coords(coords)

        _, r_ij, indices, _ = self._simulation_box._lattice.get_points_in_sphere(
                self._frac_coords, cartesian_coords,
                self.cutoff, zip_results=False,
                )

        # Need to check for length of r_ij because of
        # https://github.com/materialsproject/pymatgen/issues/3794

        if len(r_ij) == 0:  # No neighbor, zero energy.
            return 0.

        if np.any(r_ij < 1e-3):  # Close contact, infinite energy.
            return np.inf

        es_j = self._lj_params[indices]
        x = (0.5 * (es_j[:, 1] + self.sigma)) / r_ij
        e = 4 * np.sqrt(es_j[:, 0] * self.epsilon)

        return np.sum(e * (x**12 - x**6))


def voxels_from_file(
        cif_pathname,
        grid_size=GRID_SIZE,
        *,
        cutoff=CUTOFF,
        epsilon=EPSILON,
        sigma=SIGMA,
        cubic_box=CUBIC_BOX,
        n_jobs=N_JOBS,
        only_voxels=True,
        ):
    r"""
    Calculate and return voxels from ``.cif`` file.

    Parameters
    ----------
    cif_pathname : str
       Pathname to the ``.cif`` file.
    only_voxels : bool, default=True
        Determines ``out`` type.
        
    Returns
    -------
    out : array or :class:`Grid`
        If ``only_voxels=True`` array, else :class:`Grid`.

    See Also
    --------
    :func:`voxels_from_dir`
        For a description of the parameters.
    """
    grid = Grid(grid_size, cutoff=cutoff, epsilon=epsilon, sigma=sigma)

    grid.load_structure(cif_pathname)
    grid.calculate(cubic_box=cubic_box, n_jobs=n_jobs)

    if only_voxels:
        return grid.voxels

    return grid


def voxels_from_files(
        cif_pathnames,
        out_pathname,
        grid_size=GRID_SIZE,
        *,
        cutoff=CUTOFF,
        epsilon=EPSILON,
        sigma=SIGMA,
        cubic_box=CUBIC_BOX,
        n_jobs=N_JOBS,
        ):
    r"""
    Calculate voxels from a list of ``.cif`` files and store them.

    Voxels are stored under ``out_pathname`` as ``.npy`` files.

    Parameters
    ----------
    cif_pathnames : list
       List of pathnames to the ``.cif`` files.
    out_pathname : str
        Pathname to the directory under which voxels are stored.

    See Also
    --------
    :func:`voxels_from_dir`
        For a description of the parameters.

    Notes
    -----
    Structures that can't be processsed are omitted.
    """
    os.mkdir(out_pathname)

    for file in tqdm(cif_pathnames, desc='Creating energy voxels'):
        try:
            name = Path(file).stem  # Name of the structure.
            grid = voxels_from_file(
                    cif_pathname=file,
                    grid_size=grid_size,
                    cutoff=cutoff,
                    epsilon=epsilon,
                    sigma=sigma,
                    cubic_box=cubic_box,
                    n_jobs=n_jobs,
                    )

            pathname = os.path.join(out_pathname, name)
            np.save(pathname, grid)  # .npy extension is appended by default.
        except Exception as e:
            print(e)


def voxels_from_dir(
        cif_dirname: str,
        out_pathname: str,
        grid_size: int = GRID_SIZE,
        *,
        cutoff: float = CUTOFF,
        epsilon: float = EPSILON,
        sigma: float = SIGMA,
        cubic_box: float | None = CUBIC_BOX,
        n_jobs: int | None = N_JOBS,
        ):
    r"""
    Calculate voxels from a directory of ``.cif`` files and store them.

    Voxels are stored under ``out_pathname`` as ``.npy`` files.

    Parameters
    ----------
    cif_dirname : str
       Pathname to the directory containing the ``.cif`` files.
    out_pathname : str
        Pathname to the directory under which voxels are stored.
    grid_size : int, default=25
        Number of grid points along each dimension.
    cutoff : float, default=10.0
        Cutoff radius (Å) for the LJ potential.
    epsilon : float, default=50.0
        Epsilon value (ε/K) of the probe atom.
    sigma : float, default=2.5
        Sigma value (σ/Å) of the probe atom.
    cubic_box : float or None, default=None
        If ``None``, the simulation box is a supercell scaled according to
        MIC. Otherwise, cubic box of size ``cubic_box``.
    n_jobs : int, optional
        Number of jobs to run in parallel. If ``None``, then the number returned
        by ``os.cpu_count()`` is used.

    Notes
    -----
    Structures that can't be processsed are omitted.
    """
    cif_pathanmes = [os.path.join(cif_dirname, f) for f in os.listdir(cif_dirname)]

    voxels_from_files(
            cif_pathanmes, out_pathname,
            grid_size=grid_size,
            cutoff=cutoff,
            epsilon=epsilon,
            sigma=sigma,
            cubic_box=cubic_box,
            n_jobs=n_jobs,
            )
