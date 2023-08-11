import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def plot_voxels_mpl(voxels, *, fill_pattern=None, colorbar=True, cmap='viridis', **kwargs):
    r"""
    Visualizing voxels with `Axes3d.voxels
    <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.voxels.html#mpl_toolkits.mplot3d.axes3d.Axes3D.voxels>`_.

    Parameters
    ----------
    voxels : 3D array
    fill_pattern : 3D array of bool, optional
        A 3D array of values, with truthy values indicating which voxels to
        fill.
    colorbar : bool, default=True
        Whether to include a colorbar.
    cmap : str, default='viridis'
        `Colormap
        <https://matplotlib.org/stable/tutorials/colors/colormaps.html>`_ that
        colorizes the voxels based on their value.
    **kwargs :
        Valid keyword arguments for `Axes3d.voxels
        <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.voxels.html#mpl_toolkits.mplot3d.axes3d.Axes3D.voxels>`_.

        .. warning::
            Do not pass the argument ``facecolors`` of `Axes3d.voxels
            <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.voxels.html#mpl_toolkits.mplot3d.axes3d.Axes3D.voxels>`_.
            This argument is used under the hood by :func:`plot_voxels_mpl` to
            generate the colors of the voxels based on the specified ``cmap``.

    Returns
    -------
    fig : `Figure <https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure>`_
    """
    if fill_pattern == None:
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


def plot_voxels_pv():
    pass
