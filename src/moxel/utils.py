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
import json
import itertools
import numpy as np
from pathlib import Path
import matplotlib as mpl
from itertools import repeat
from rich.progress import track
import matplotlib.pyplot as plt
from multiprocessing import Pool
from dataclasses import dataclass
from pymatgen.core import Structure


def mic_scale_factors(r, lattice_vectors):
    r"""
    Return scale factors to satisfy minimum image convention (MIC) [1]_ .

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
    .. [1] W. Smith, "The Minimum Image Convention in Non-Cubic MD Cells", 1989.
    """
    a, b, c = lattice_vectors
    volume = np.linalg.norm(np.dot(a, np.cross(b, c)))

    w_a = volume/np.linalg.norm(np.cross(b, c))
    w_b = volume/np.linalg.norm(np.cross(a, c))
    w_c = volume/np.linalg.norm(np.cross(a, b))

    return np.ceil(2*r/np.array([w_a, w_b, w_c]))


#def mic_scale_factors(r, lattice_vectors):
#    r"""
#    Return scale factors to satisfy minimum image convention (MIC).
#
#    Parameters
#    ----------
#    r : float
#        The cutoff radius used in MIC convetion.
#    lattice_vectors : array of shape (3, 3)
#        The lattice vectors of the unit cell.
#        Each row corresponds to a lattice vector.
#
#    Returns
#    -------
#    scale_factors : array of shape (3,)
#        ``scale_factors[i]`` scales ``lattice_vectors[i]``.
#    """
#    a, b, c = lattice_vectors
#
#    a_hat = a/np.linalg.norm(a)
#    b_a = b - np.dot(b, a_hat)*a_hat # Rejection of b from a.
#
#    ab_cross = np.cross(a, b)
#    ab_cross_hat = ab_cross/np.linalg.norm(ab_cross)
#    c_ab = np.dot(c, ab_cross_hat)
#
#    norms = [np.linalg.norm(i) for i in [a, b_a, c_ab]]
#
#    return np.ceil(2*r/np.array(norms))


@dataclass
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
    structure : `pymatgen.core.structure.Structure <https://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.Structure>`_
        Available only after :meth:`Grid.load_structure` has been called.
    structure_name : str
        Available only after :meth:`Grid.load_structure` has been called.
    cubic_box : bool
        Available only after :meth:`Grid.calculate` has been called.
    voxels : array of shape (grid_size, grid_size, grid_size)
       Available only after :meth:`Grid.calculate` has been called.
    """
    def __init__(self, grid_size=25, cutoff=10, epsilon=50, sigma=2.5):
        self.grid_size = grid_size
        self.cutoff = cutoff
        self.epsilon = epsilon
        self.sigma = sigma

    def load_structure(self, cif_pathname):
        r"""
        Load a crystal structure from a ``.cif`` file.

        Parameters
        ----------
        cif_pathname : str
           Pathname to the ``.cif`` file.
        """
        self.structure = Structure.from_file(cif_pathname)
        #self.structure_name = cif_pathname.split('/')[-1].split('.')[0]
        self.structure_name = cif_pathname.split('/')[-1].removesuffix('.cif')

    def calculate(self, cubic_box=False, length=30, potential='lj'):
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
            The size of the cubic box in Å. Takes effect only if ``cubic_box=True``.

        Returns
        -------
        voxels : array of shape (grid_size, grid_size, grid_size)
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
            probe_coords = np.linspace(0-d, 0+d, self.grid_size) # Cartesian
            self._simulation_box = self.structure
        else:
            probe_coords = np.linspace(0, 1, self.grid_size) # Fractional
            scale = mic_scale_factors(self.cutoff, self.structure.lattice.matrix)
            self._simulation_box = self.structure * scale
        
        if potential == 'lj':
            # Embarrassingly parallel.
            with Pool(processes=None) as p:
                energies = list(
                        p.map(
                            self.lj_potential, itertools.product(*(probe_coords,)*3)
                            )
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
            If ``cubic_box=True`` cartesian. Else, fractional.

        Returns
        -------
        energy : float
            Energy as :math:`e^{-\beta \mathcal{V}}`, to ensure numerical stability.
        """
        if self.cubic_box:
            cartesian_coords = coords
            neighbors = self._simulation_box.get_sites_in_sphere(coords, self.cutoff)
        else:
            cartesian_coords = self._simulation_box.lattice.get_cartesian_coords(coords)
            neighbors = self._simulation_box.get_sites_in_sphere(cartesian_coords, self.cutoff)

        energy = 0
        if len(neighbors) != 0:
            for atom in neighbors:
                r_ij = np.linalg.norm(cartesian_coords - atom.coords)
                if r_ij <= 1e-3:
                    energy += 1000
                else:
                    e_j, s_j = lj_params[atom.species_string]
                    x = (0.5 * (s_j + self.sigma)) / r_ij
                    energy += 4 * np.sqrt(e_j * self.epsilon) * (x**12 - x**6)

        return np.exp(-(1 / 298) * energy) # For numerical stability.


def voxels_from_file(
        cif_pathname, grid_size=25, cutoff=10,
        epsilon=50, sigma=2.5, cubic_box=False, length=30, 
        only_voxels=True,
        ):
    """
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
        The size of the cubic box in Å. Takes effect only if ``cubic_box=True``.
    only_voxels : bool, default=True
        Determines ``out`` type.
        
    Returns
    -------
    out : ``array`` or :class:`.Grid`
        If ``only_voxels=True``, array of shape ``(grid_size, grid_size,
        grid_size)``. Otherwise, :class:`.Grid`.

    Notes
    -----

    * For structures that can not be processsed, their voxels are filled with zeros.
    """
    grid = Grid(grid_size, cutoff, epsilon, sigma)
    try:
        grid.load_structure(cif_pathname)
        grid.calculate(cubic_box=cubic_box, length=length)
    except ValueError:
        grid.voxels = np.full(shape=(grid_size,)*3, fill_value=0, dtype=np.float32)

    if only_voxels:
        return grid.voxels
    else:
        return grid


def voxels_from_files(
        cif_pathnames, grid_size=25, cutoff=10,
        epsilon=50, sigma=2.5,
        cubic_box=False, length=30,
        out_pathname=None,
        ):
    """
    Calculate voxels from a list of ``.cif`` files and save them in ``out_pathname`` as
    ``array`` of shape ``(n_samples, grid_size, grid_size, grid_size)``, where
    ``n_samples`` is the number of is the number of ``.cif`` files in
    ``cif_pathnames``.


    Parameters
    ----------
    cif_pathanmes : list
       List of pathnames to the ``.cif`` files.
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
        The size of the cubic box in Å. Takes effect only if ``cubic_box=True``.
    out_pathname : str, optional
        Pathname to the file holding the voxels. If not specified, voxels are stored in
        ``./voxels.npy``.
    
    Notes
    -----

    * Samples in output array follow the order in ``cif_pathnames``.

    * For structures that can not be processsed, their voxels are filled with zeros.
    """
    n = len(cif_pathnames)
    cif_files = [i for i in cif_pathnames if i.endswith('.cif')]
    out_pathname = './voxels.npy' if not out_pathname else out_pathname

    fp = np.lib.format.open_memmap(
        out_pathname, mode='w+',
        shape=(n, *(grid_size,)*3),
        dtype=np.float32
        )

    grids = map(
            voxels_from_file, cif_files,
            repeat(grid_size), repeat(cutoff),
            repeat(epsilon), repeat(sigma),
            repeat(cubic_box), repeat(length)
            )

    for i in track(range(n), description='Processing...'):
        fp[i] = next(grids)

    fp.flush()


def voxels_from_dir(
        cif_dirname, grid_size=25, cutoff=10,
        epsilon=50, sigma=2.5,
        cubic_box=False, length=30,
        out_pathname=None,
        ):
    """
    Calculate voxels from ``.cif`` files in ``cif_dirname`` and save them in
    ``out_pathname`` as ``array`` of shape ``(n_samples, grid_size, grid_size,
    grid_size)``, where ``n_samples`` is the number of ``.cif`` files in
    ``cif_dirname``.

    Parameters
    ----------
    cif_dirname : str
       Pathname to the directory containing the ``.cif`` files.
    grid_size : int, default=25
       Number of grid points along each dimension.
    cutoff : float, default=10
        Cutoff radius (Å) for the LJ potential.
    epsilon : float, default=50
        Epsilon value (ε/K) of the probe atom.
    cubic_box : bool, default=False
        If ``True``, the simulation box is cubic.
    length : float, default=30
        The size of the cubic box in Å. Takes effect only if ``cubic_box=True``.
    sigma : float, default=25
        Sigma value (σ/Å) of the probe atom.
    out_pathname : str, optional
        Pathname to the file holding the voxels. If not specified, voxels are stored
        in ``./voxels.npy``.
    
    Notes
    -----

    * Samples in output array follow the order in ``sorted(os.listdir(cif_dirname))``.

    * For structures that can not be processsed, their voxels are filled with zeros.
    """
    files = sorted(os.listdir(cif_dirname))
    cif_files = [f'{cif_dirname}/{file}' for file in files if file.endswith('.cif')]
    n = len(cif_files)
    out_pathname = './voxels.npy' if not out_pathname else out_pathname

    fp = np.lib.format.open_memmap(
            out_pathname, mode='w+',
            shape=(n, *(grid_size,)*3),
            dtype=np.float32
            )

    grids = map(
            voxels_from_file, cif_files,
            repeat(grid_size), repeat(cutoff),
            repeat(epsilon), repeat(sigma),
            repeat(cubic_box), repeat(length)
            )

    for i in track(range(n), description='Processing...'):
        fp[i] = next(grids)

    fp.flush()


#def batch_create(cif_pathnames, n_batches, batches_dirname=None):
#    """
#    Split a number of structures into ``n_batches`` of approximately equal size.
#
#    The batches are created under the ``batches_dirname`` directory.
#
#    For example, the i-th batch corresponds to ``batches_dirname/batch_i``.
#
#    Note that the **structures are randomly shuffled** prior to splitting.  The
#    new ordering for the i-th batch is stored in
#    ``batches_dirname/batch_i/names.json``
#
#    Parameters
#    ----------
#    cif_filenames : list
#        List of pathnames to the ``.cif`` files.
#    n_batches : int
#        Number of batches that will be created.
#    batches_dirname : str, optional
#        Pathname to the directory where batches will be created. If ``None``, the
#        batches are created under ``./``.
#   """
#    cif_pathnames = np.array([i for i in cif_pathnames if i.endswith('.cif')])
#    np.random.shuffle(cif_pathnames)
#
#    batches = np.array_split(cif_pathnames, n_batches)
#
#    if batches_dirname == None:
#        batches_dirname = '.'
#
#    for i, batch in enumerate(batches):
#        os.mkdir(f'{batches_dirname}/batch_{i}')
#        info_dict = {'names': list(batch), 'size': len(batch)}
#
#        with open(f'{batches_dirname}/batch_{i}/names.json', 'w') as fhand:
#            json.dump(info_dict, fhand, indent=4)
#

#def batch_calculate(batch_dirname, cif_dirname, **kwargs):
#    """
#    Calculate voxels for the structures in ``batch_dirname/names.json``.
#
#    The voxels are saved in ``batch_dirname/voxels.npy``.
#
#    Parameters
#    ----------
#    batch_dirname : str
#        Pathname to the batch directory.
#    cif_dirname : str
#        Pathname to the directory holding the ``.cif`` files.
#    **kwargs :
#        Valid keyword arguments for :func:`voxels_from_files`.
#
#        .. warning::
#            Do not pass the arguments ``cif_pathnames`` and ``out_pathname`` of
#            :func:`voxels_from_files`.
#
#    Notes
#    -----
#    For structures that can not be processsed, their voxels are filled with
#    zeros.
#    """
#    with open(f'{batch_dirname}/names.json', 'r') as fhand:
#        info = json.load(fhand)
#
#    names, size = info['names'], info['size']
#    cif_names = [f'{cif_dirname}/{name}.cif' for name in names]
#
#    voxels_from_files(
#        cif_names,
#        out_pathname=f'{batch_dirname}/voxels.npy',
#        **kwargs
#        )
#

def batch_clean_and_merge(batch_dirnames, out_pathname=None):
    """
    Clean a single batch, or *first clean and then merge* multiple batches.
    All batches must have the form::

        batch
        ├──voxels.npy
        └──names.json

    Cleaning is required since the voxels for some structures might be zero,
    see :func:`Grid.calculate`.

    If ``len(batch_dirnames) == 1`` the cleaned voxels for are stored under
    ``batch_dirnames[0]/clean_voxels.npy`` and the names of their corresponding
    structures under ``batch_dirnames[0]/clean_names.json``.

    If ``len(batch_dirnames) > 1`` the voxels (*cleaned and merged*) are stored
    under ``out_pathname/clean_voxels.npy`` and the names of their corresponding
    structures under ``out_pathname/clean_names.json``. That is::

        out_pathname
        ├──clean_voxels.npy
        └──clean_names.json

    Parameters
    ----------
    batch_dirnames : list
        List of pathnames to the directories of the batches.
    out_pathname : str, optional
        Pathname to the directory holding the clean voxels and names.  The
        directory is created if it doesn't exist. Takes effect only if
        ``len(batch_dirnames) > 1``. If ``None`` voxels and names are
        stored under ``./clean_voxels.npy`` and ``./clean_names.json``.

    Returns
    -------
    exit_status : int
        If no voxels are missing ``0`` else ``1``.
    """
    batch_dict = dict()

    for i, batch_dir in enumerate(batch_dirnames):
        with open(f'{batch_dir}/names.json', 'r') as fhand:
            info = json.load(fhand)

        names = np.array(info['names'])

        fp = np.load(f'{batch_dir}/voxels.npy', mmap_mode='r')

        missing_idx = [i for i, j in enumerate(fp) if np.all(j == 0)]

        batch_dict[f'batch_{i}'] = (fp, names, missing_idx)

    missing = sum([len(batch_dict[i][2]) != 0 for i in batch_dict])
    if missing == 0:
        os.rename(f'{batch_dirnames[0]}/names.json', f'{batch_dirnames[0]}/clean_names.json') 
        os.rename(f'{batch_dirnames[0]}/voxels.npy', f'{batch_dirnames[0]}/clean_voxels.npy') 
        print('No missing voxels found!!')
        return 0

    print('Missing voxels found! Cleaning...')
    if len(batch_dirnames) == 1:
        clean_dir = batch_dirnames[0]
    elif out_pathname == None:
        clean_dir = '.'
    else:
        clean_dir = out_pathname
        try:
            os.mkdir(clean_dir)
        except:
            pass

    clean_size = np.sum([len(batch_dict[b][0]) - len(batch_dict[b][2]) for b in batch_dict.keys()])
    clean_fp = np.lib.format.open_memmap(
        f'{clean_dir}/clean_voxels.npy',
        shape=(clean_size, *fp.shape[1:]),
        mode='w+', dtype='float32',
        )

    clean_idx = 0
    for b in batch_dict.keys():
        for idx, x in enumerate(batch_dict[b][0]):
            if idx in batch_dict[b][2]:
                pass
            else:
                clean_fp[clean_idx] = x
                clean_idx += 1

    clean_names = np.concatenate([np.delete(batch_dict[b][1], batch_dict[b][2]) for b in batch_dict.keys()])
    clean_dict = {'names': list(clean_names)}

    with open(f'{clean_dir}/clean_names.json', 'w') as fhand:
        json.dump(clean_dict, fhand, indent=4)

    return 1


def plot_voxels(voxels, *, fill_pattern=None, colorbar=True, cmap='viridis', **kwargs):
    r"""
    Visualizing voxels with `Axes3d.voxels`_.

    .. _Axes3d.voxels: https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.voxels.html#mpl_toolkits.mplot3d.axes3d.Axes3D.voxels

    Parameters
    ----------
    voxels : 3D array
    fill_pattern : 3D array of bool, optional
        A 3D array of truthy values, indicating which voxels to fill. If not
        specified, all voxels are filled.
    colorbar : bool, default=True
        Whether to include a colorbar.
    cmap : str, default='viridis'
        `Colormap
        <https://matplotlib.org/stable/tutorials/colors/colormaps.html>`_ that
        colorizes the voxels based on their value.
    **kwargs :
        Valid keyword arguments for `Axes3d.voxels`_.

        .. warning::
            Do not pass the argument ``facecolors`` of `Axes3d.voxels`_.
            This argument is used under the hood by :func:`plot_voxels_mpl` to
            generate the colors of the voxels based on the specified ``cmap``.

    Returns
    -------
    fig : `Figure <https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure>`_
    """

    if np.all(fill_pattern) == None:
        fill_pattern = np.full(voxels.shape, True)

    cmap = plt.get_cmap(cmap)
    norm = mpl.colors.Normalize()
    colors = cmap((voxels - voxels.min()) / (voxels.max() - voxels.min()))

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    _ = ax.voxels(filled=fill_pattern, facecolors=colors, **kwargs)

    if colorbar:
        mappable = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(), cmap=cmap) 
        cbar = fig.colorbar(
            mappable, ax=ax, ticks=[], extend='max',
            shrink=0.7, pad=0.1,
            )
        cbar.ax.set_ylabel(r'$e^{-\beta \mathcal{V}}$', fontsize=12)
    
    return fig


# Loading LJ parameters.
with open(f'{Path(__file__).parents[0]}/lj_params.json', 'r') as fhand:
    lj_params = json.load(fhand)
