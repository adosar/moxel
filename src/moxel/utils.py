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
"""

import os
import json
import itertools
from pathlib import Path
from multiprocessing import Pool
from itertools import repeat
import warnings
import numpy as np
from scipy.special import erfc
from tqdm import tqdm
from pymatgen.core import Structure
from . _params import lj_params
warnings.filterwarnings('ignore')


def get_names(fname):
    r"""
    Load a list of material names saved in ``.json`` format.

    Parameters
    ----------
    fname : str
        Pathname to the ``.json`` file.

    Returns
    -------
    names : list
    """
    with open(fname, 'r') as fhand:
        names = json.load(fhand)

    return names


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

    return np.ceil(2*r/np.array([w_a, w_b, w_c]))


class Grid:
    r"""
    A 3D energy grid over a crystal structure.

    Parameters
    ----------
    grid_size : int, default=25
        Number of grid points along each dimension.
    cutoff : float, default=10
        Cutoff radius (Å) for the LJ potential.
    epsilon : float, default=50
        Epsilon value (ε/K) of the probe atom.
    sigma : float, default=2.5
        Sigma value (σ/Å) of the probe atom.

    Attributes
    ----------
    structure : :class:`pymatgen.core.structure.Structure`
        Available only after :meth:`Grid.load_structure` has been called.
    structure_name : str
        Available only after :meth:`Grid.load_structure` has been called.
    cubic_box : bool
        Available only after :meth:`Grid.calculate` has been called.
    voxels : array of shape (grid_size,)*3
       Available only after :meth:`Grid.calculate` has been called.
    """
    def __init__(self, grid_size=25, cutoff=10, epsilon=50, sigma=2.5, wolf_alpha=0.2, charge_prop="charge"):
        self.grid_size = grid_size
        self.cutoff = cutoff
        self.epsilon = epsilon
        self.sigma = sigma
        self.wolf_alpha = wolf_alpha
        self.charge_prop = charge_prop

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

    def calculate(self, cubic_box=False, length=30, potential='lj', n_jobs=None):
        r"""
        Iterate over the grid and return voxels.

        For computational efficiency and to assure (approximately) the same
        spatial resolution, the grid is overlayed over a supercell scaled
        according to MIC, see :func:`mic_scale_factors`.

        If lattice angles are significantly different than 90°, to avoid
        distortions set ``cubic_box`` to ``True``. In this case, the grid is
        overlayed over a cubic box of size ``length`` centered at the origin but
        periodicity is no longer guaranteed.

        Parameters
        ----------
        potential : str, default='lj'
            The potential used to calculate voxels. Currently, only the
            LJ potential is supported.
        cubic_box : bool, default=False
            If ``True``, the simulation box is cubic.
        length : float, default=30
            The size of the cubic box in Å. Takes effect only
            if ``cubic_box == True``.
        n_jobs : int, optional
            Number of jobs to run in parallel. If ``None``, all processors are
            used.

        Returns
        -------
        voxels : array of shape (grid_size,)*3
            The energy voxels as :math:`e^{-\beta \mathcal{V}}`, to ensure
            numerical stability.
        Notes
        -----
        For structures that can not be processsed, their voxels are filled with
        zeros.
        """
        self.cubic_box = cubic_box

        if cubic_box:
            d = length / 2
            probe_coords = np.linspace(0-d, 0+d, self.grid_size, endpoint=False)  # Cartesian.
            self._simulation_box = self.structure
        else:
            probe_coords = np.linspace(0, 1, self.grid_size, endpoint=False)  # Fractional.
            scale = mic_scale_factors(self.cutoff, self.structure.lattice.matrix)
            self._simulation_box = self.structure * scale

        func = None
        if potential == 'lj':
            func = self.lj_potential
            # Cache LJ parameters for all atoms in the simulation box.
            self._lj_params = np.array([lj_params[atom.species_string] for atom in self._simulation_box])
        else:
            if potential == 'elecpot':
                func = self.elec_potential
            elif potential == 'elecfield':
                func = self.elec_field
            elif potential == 'elecfield_vect':
                func = self.elec_field_vect
            try:
                # Cache charges for all atoms in the simulation box
                self._charges = np.array(self._simulation_box.site_properties[self.charge_prop])
            except KeyError as e:
                raise ValueError("atomic charges could not be found") from e

        if func is None:
            raise ValueError(f"invalid potential type: '{potential}'")

        # Cache fractional coordinates, since this is a slow function in pymatgen.
        self._frac_coords = self._simulation_box.frac_coords
        # Embarrassingly parallel.
        with Pool(processes=n_jobs) as p:
            energies = p.map(func, itertools.product(*(probe_coords,)*3))

        if potential == 'elecfield_vect':
            self.voxels = np.array(energies, dtype=np.float32).reshape((self.grid_size, self.grid_size, self.grid_size, 3))
        else:
            self.voxels = np.array(energies, dtype=np.float32).reshape((self.grid_size,)*3)

        return self.voxels

    def lj_potential(self, coords):
        r"""
        Calculate LJ potential at cartesian or fractional
        coordinates.

        Parameters
        ----------
        coordinates : array_like of shape (3,)
            If ``cubic_box == True`` cartesian. Else, fractional.

        Returns
        -------
        energy : float
            Energy as :math:`e^{-\beta \mathcal{V}}`, to ensure numerical stability.
        """
        if self.cubic_box:
            cartesian_coords = coords
        else:
            cartesian_coords = self._simulation_box.lattice.get_cartesian_coords(coords)

        _, r_ij, indices, _ = self._simulation_box._lattice.get_points_in_sphere(
                self._frac_coords, cartesian_coords,
                self.cutoff, zip_results=False,
                )

        '''
        Need to check for length of r_ij because of
        https://github.com/materialsproject/pymatgen/issues/3794
        '''
        if len(r_ij) == 0:  # No neighbor, zero energy.
            return 1.

        if np.any(r_ij < 1e-3):  # Close contact, infinite energy.
            return 0.

        es_j = self._lj_params[indices]
        x = (0.5 * (es_j[:, 1] + self.sigma)) / r_ij
        e = 4 * np.sqrt(es_j[:, 0] * self.epsilon)
        energy = sum(e * (x**12 - x**6))

        # This should be changed with clipping in future versions.
        return np.exp(-(1 / 298) * energy)  # For numerical stability.

    def elec_potential(self, coords):
        r"""
        Calculate electrostatic potential at cartesian or fractional
        coordinates, in the Wolf summation approximation, in its
        damped shifted potential (DSP) formulation. For details, see
        the PhD thesis from Guillaume Fraux, section 6.3.3, at
        https://theses.fr/2019PSLEC004

        Parameters
        ----------
        coordinates : array_like of shape (3,)
            If ``cubic_box == True`` cartesian. Else, fractional.

        Returns
        -------
        potential: float
            Electric potential.
        """
        if self.cubic_box:
            cartesian_coords = coords
        else:
            cartesian_coords = self._simulation_box.lattice.get_cartesian_coords(coords)

        _, r_ij, indices, _ = self._simulation_box._lattice.get_points_in_sphere(self._frac_coords, cartesian_coords, self.cutoff, zip_results=False)

        # No neighbor, zero potential
        # Need to check for length of r_ij because of https://github.com/materialsproject/pymatgen/issues/3794
        if len(r_ij) == 0:
            return 0.

        # Implement Eq. (6.31) from https://theses.fr/2019PSLEC004
        q_j = self._charges[indices]
        sum_q = np.sum(q_j)
        pot = q_j * erfc(self.wolf_alpha * r_ij) / r_ij
        z = erfc(self.wolf_alpha * self.cutoff) / self.cutoff
        return np.sum(pot) - sum_q * z + (z / 2 + self.wolf_alpha / np.sqrt(np.pi))

    def elec_field_vect(self, coords):
        r"""
        Calculate electrostatic field at cartesian or fractional
        coordinates, in the Wolf summation approximation, in its
        damped shifted force (DSF) formulation. For details, see
        the PhD thesis from Guillaume Fraux, section 6.3.3, at
        https://theses.fr/2019PSLEC004

        Parameters
        ----------
        coordinates : array_like of shape (3,)
            If ``cubic_box == True`` cartesian. Else, fractional.

        Returns
        -------
        field: numpy.ndarray of shape (3,)
            Electric field vector.
        """
        if self.cubic_box:
            cartesian_coords = coords
        else:
            cartesian_coords = self._simulation_box.lattice.get_cartesian_coords(coords)

        frac_j, r_ij, indices, _ = self._simulation_box._lattice.get_points_in_sphere(self._frac_coords, cartesian_coords, self.cutoff, zip_results=False)

        # No neighbor, zero field
        # Need to check for length of r_ij because of https://github.com/materialsproject/pymatgen/issues/3794
        if len(r_ij) == 0:
            return 0.

        # We implement Eq. (6.35) from https://theses.fr/2019PSLEC004
        # fac = the term the term inside the brackets
        q_j = self._charges[indices]
        r2_ij = r_ij**2
        fac = erfc(self.wolf_alpha * r_ij) / r2_ij + 2 * self.wolf_alpha/np.sqrt(np.pi) * np.exp(-self.wolf_alpha**2 * r2_ij) / r_ij
        fac -= erfc(self.wolf_alpha * self.cutoff) / self.cutoff**2 + 2 * self.wolf_alpha/np.sqrt(np.pi) * np.exp(-self.wolf_alpha**2 * self.cutoff**2) / self.cutoff

        # Calculate the r_ij vectors, and normalize them.
        frac_j -= coords
        vec_ij = self._simulation_box._lattice.get_cartesian_coords(frac_j)
        vec_ij /= r_ij[:, np.newaxis]

        # Make sure the round-tripping of PBC is consistent between
        # get_points_in_sphere() and the handwritten code above.
        assert np.max(np.fabs(np.sum(vec_ij**2, axis=1) - 1)) < 1.e-12

        # Put everthing together.
        field = np.sum(q_j[:, np.newaxis] * fac[:, np.newaxis] * vec_ij, axis=0)
        return field

    def elec_field(self, coords):
        r"""
        Calculate electrostatic field strength at cartesian or fractional
        coordinates, in the Wolf summation approximation. Same as
        ``elec_field_vect``, but returns the norm of the field as a
        scalar.

        Parameters
        ----------
        coordinates : array_like of shape (3,)
            If ``cubic_box == True`` cartesian. Else, fractional.

        Returns
        -------
        field: scalar
            Electric field norm.
        """
        field = self.elec_field_vect(coords)
        return np.sqrt(np.sum(field**2))


def voxels_from_file(
        cif_pathname, grid_size=25, cutoff=10,
        epsilon=50, sigma=2.5, cubic_box=False, length=30, 
        n_jobs=None, only_voxels=True,
        ):
    r"""
    Return voxels from ``.cif`` file.

    Parameters
    ----------
    cif_pathname : str
       Pathname to the ``.cif`` file.
    grid_size : int, default=25
        Number of grid points along each dimension.
    cutoff : float, default=10
        Cutoff radius (Å) for the LJ potential.
    epsilon : float, default=50
        Epsilon value (ε/K) of the probe atom.
    sigma : float, default=25
        Sigma value (σ/Å) of the probe atom.
    cubic_box : bool, default=False
        If ``True``, the simulation box is cubic.
    length : float, default=30
        The size of the cubic box in Å. Takes effect only if ``cubic_box == True``.
    n_jobs : int, optional
        Number of jobs to run in parallel. If ``None``, all processors are used.
    only_voxels : bool, default=True
        Determines ``out`` type.
        
    Returns
    -------
    out : ``array`` or :class:`Grid`
        If ``only_voxels == True``, array of shape ``(grid_size,)*3``.
        Otherwise, :class:`Grid`.

    Notes
    -----
    For structures that can not be processsed, their voxels are filled with zeros.
    """
    grid = Grid(grid_size, cutoff, epsilon, sigma, charge_prop="pbe_ddec_charge")
    try:
        grid.load_structure(cif_pathname)
        grid.calculate(cubic_box=cubic_box, length=length, n_jobs=n_jobs, potential="elecpot")
    except:
        grid.voxels = np.full(shape=(grid_size,)*3, fill_value=0, dtype=np.float32)

    if only_voxels:
        return grid.voxels
    else:
        return grid


def voxels_from_files(
        cif_pathnames, out_pathname, grid_size=25, cutoff=10,
        epsilon=50, sigma=2.5, cubic_box=False, length=30,
        n_jobs=None,
        ):
    r"""
    Calculate voxels from a list of ``.cif`` files and store them under
    ``out_pathname`` as :class:`numpy.array` of shape
    ``(n_samples, grid_size, grid_size, grid_size)``,
    where ``n_samples == len(cif_pathnames)``.

    After processing the following files are created::

        out_pathname
            ├──voxels.npy
            └──names.json

    The file ``names.json`` stores the names of the materials as a
    :class:`list`, which might be useful for later indexing.

    Parameters
    ----------
    cif_pathnames : list
       List of pathnames to the ``.cif`` files.
    out_pathname : str
        Pathname to the directory under which voxels are stored.
    grid_size : int, default=25
        Number of grid points along each dimension.
    cutoff : float, default=10
        Cutoff radius (Å) for the LJ potential.
    epsilon : float, default=50
        Epsilon value (ε/K) of the probe atom.
    sigma : float, default=25
        Sigma value (σ/Å) of the probe atom.
    cubic_box : bool, default=False
        If ``True``, the simulation box is cubic.
    length : float, default=30
        The size of the cubic box in Å. Takes effect only if ``cubic_box == True``.
    n_jobs : int, optional
        Number of jobs to run in parallel. If ``None``, then the number returned
        by ``os.cpu_count()`` is used.

    Notes
    -----
    * Samples in output array follow the order in ``cif_pathnames``.
    * For structures that can not be processsed, their voxels are filled with zeros.
    """
    n_samples = len(cif_pathnames)
    names = [Path(i).stem for i in cif_pathnames]

    # Store the names.
    with open(f'{out_pathname}/names.json', mode='w') as fhand:
        json.dump(names, fhand, indent=4)

    fp = np.lib.format.open_memmap(
        f'{out_pathname}/voxels.npy', mode='w+',
        shape=(n_samples, *(grid_size,)*3),
        dtype=np.float32
        )

    grids = map(
            voxels_from_file, cif_pathnames,
            repeat(grid_size), repeat(cutoff),
            repeat(epsilon), repeat(sigma),
            repeat(cubic_box), repeat(length),
            repeat(n_jobs)
            )

    for i in tqdm(range(n_samples), desc='Creating voxels'):
        fp[i] = next(grids)

    fp.flush()


def voxels_from_dir(
        cif_dirname, out_pathname, grid_size=25, cutoff=10,
        epsilon=50, sigma=2.5, cubic_box=False, length=30,
        n_jobs=None,
        ):
    r"""
    Calculate voxels from a directory of ``.cif`` files and save them under
    ``out_pathname`` as :class:`numpy.array` of shape
    ``(n_samples, grid_size, grid_size, grid_size)``,
    where ``n_samples == len(cif_pathnames)``.

    After processing the following files are created::

        out_pathname
            ├──voxels.npy
            └──names.json

    The file ``names.json`` stores the names of the materials as a
    :class:`list`, which might be useful for later indexing.

    Parameters
    ----------
    cif_dirname : str
       Pathname to the directory containing the ``.cif`` files.
    out_pathname : str
        Pathname to the directory under which voxels are stored.
    grid_size : int, default=25
        Number of grid points along each dimension.
    cutoff : float, default=10
        Cutoff radius (Å) for the LJ potential.
    epsilon : float, default=50
        Epsilon value (ε/K) of the probe atom.
    sigma : float, default=25
        Sigma value (σ/Å) of the probe atom.
    cubic_box : bool, default=False
        If ``True``, the simulation box is cubic.
    length : float, default=30
        The size of the cubic box in Å. Takes effect only if ``cubic_box == True``.
    n_jobs : int, optional
        Number of jobs to run in parallel. If ``None``, then the number returned
        by ``os.cpu_count()`` is used.

    Notes
    -----
    * Samples in output array follow the order in ``sorted(os.listdir(cif_dirname))``.
    * For structures that can not be processsed, their voxels are filled with zeros.
    """
    files = [os.path.join(cif_dirname, f) for f in sorted(os.listdir(cif_dirname))]

    voxels_from_files(
            files, out_pathname,
            grid_size=grid_size, cutoff=cutoff,
            epsilon=epsilon, sigma=sigma,
            cubic_box=cubic_box, length=length,
            n_jobs=n_jobs,
            )


def batch_clean(batch_dirname):
    """
    Clean a single batch.

    The batch must have the form::

        batch
        ├──voxels.npy
        └──names.json

    Cleaning is required since the voxels for some structures might be zero,
    see :func:`Grid.calculate`. After cleaning, the directory has the form::

        batch
        ├──voxels.npy
        ├──names.json
        ├──clean_voxels.npy
        └──clean_names.json

    Parameters
    ----------
    batch_dirname : str
        Pathname to the directory which requires cleaning.

    Returns
    -------
    exit_status : int
        If no voxels are missing ``0`` else ``1``.
    """
    # Case 1: no missing voxels.
    names = get_names(f'{batch_dirname}/names.json')
    voxels = np.load(f'{batch_dirname}/voxels.npy', mmap_mode='r')

    missing_idx = [i for i, x in enumerate(voxels) if np.all(x == 0)]

    if len(missing_idx) == 0:
        print('No missing voxels found!')
        return 0

    # Case 2: missing voxels.
    print('Missing voxels found! Cleaning...')

    clean_size = len(voxels) - len(missing_idx)

    # Create a new array to store the clean voxels.
    clean_fp = np.lib.format.open_memmap(
        f'{batch_dirname}/clean_voxels.npy',
        shape=(clean_size, *voxels.shape[1:]),  # Shape (N, grid, grid, grid).
        mode='w+', dtype='float32',
        )

    clean_idx = 0
    for idx, x in enumerate(voxels):
        if idx in missing_idx:
            pass
        else:
            clean_fp[clean_idx] = x
            clean_idx += 1

    clean_names = np.delete(names, missing_idx)

    with open(f'{batch_dirname}/clean_names.json', 'w') as fhand:
        json.dump(list(clean_names), fhand, indent=4)

    return 1
