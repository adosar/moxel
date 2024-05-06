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
This module provides helper functions for visualizing voxels.

.. tip::
    For faster rendering you should prefer :func:`plot_voxels_pv`.
"""

import numpy as np
import pyvista as pv
import matplotlib as mpl
import matplotlib.pyplot as plt


def plot_voxels_mpl(voxels, *, fill_pattern=None, colorbar=True, cmap='viridis', **kwargs):
    r"""
    Visualize voxels with `Axes3d.voxels`_.

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

        .. note::
            Do not pass the argument ``facecolors`` of `Axes3d.voxels`_.
            This argument is used under the hood by :func:`plot_voxels` to
            generate the colors of the voxels based on the specified ``cmap``.

    Returns
    -------
    fig : `Figure <https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure>`_

    Examples
    --------
    >>> voxels = np.random.randn(5, 5, 5)
    >>> fig = plot_voxels_mpl(voxels, cmap='coolwarm')
    >>> plt.show()  # Not needed for Jupyter.
    """

    if np.all(fill_pattern) == None:
        fill_pattern = np.full(voxels.shape, True)

    cmap = plt.get_cmap(cmap)
    norm = mpl.colors.Normalize()
    colors = cmap((voxels - voxels.min()) / (voxels.max() - voxels.min()))

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    _ = ax.voxels(filled=fill_pattern, facecolors=colors, **kwargs)

    if colorbar:
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = fig.colorbar(
            mappable, ax=ax, ticks=[], extend='max',
            shrink=0.7, pad=0.1,
            )
        cbar.ax.set_ylabel(r'$e^{-\beta \mathcal{V}}$', fontsize=12)
    
    return fig


def plot_voxels_pv(voxels, **kwargs):
    r"""
    Visualize voxels with `Plotter.add_volume`_

    .. note::
        For interactive plots in Jupyter: ``pip install "pyvista[jupyter]"``

    .. _Plotter.add_volume: https://docs.pyvista.org/version/stable/api/plotting/_autosummary/pyvista.plotter.add_volume#pyvista-plotter-add-volume

    Parameters
    ----------
    voxels : 3D array
    **kwargs : dict, optional
        Valid keyword arguments for `Plotter.add_volume`_.

    Examples
    --------
    >>> voxels = np.random.randn(5, 5, 5)
    >>> plot_voxels_pv(voxels, cmap='coolwarm')
    """
    pl = pv.Plotter()
    pl.add_volume(voxels, **kwargs)
    pl.show()
