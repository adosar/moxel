.. highlight:: python

|:rocket:| Tutorial
===================

As stated in :ref:`advantages`, all you need is a ``.cif`` file!

If you don't have one, download the following: :download:`IRMOF-1.cif<down/IRMOF-1.cif>`.

Calculation and visualization of energy voxels
----------------------------------------------

Calculation
^^^^^^^^^^^^^^^

Note that ``path/to/`` can be an absolute or relative pathname.

1. Functional interface

    .. code-block::

        >>> from moxel import voxels_from_file # Omitting .utils also works.
        >>> voxels = voxels_from_file('path/to/IRMOF-1.cif', grid_size=25)

2. Object-oriented interface

    .. code-block::

        >>> from moxel.utils import Grid # Omitting .utils also works.
        >>> grid = Grid(grid_size=25)
        >>> grid.load_structure('path/to/IRMOF-1.cif')
        >>> grid.calculate()

    .. code-block::

        >>> np.all(voxels == grid.voxels) # A sanity check.
        True

Of course, we are interested in calculating energy voxels from multiple files.
In this case, check:

* `voxels_from_files() <moxel.html#moxel.utils.voxels_from_files>`_ 
* `voxels_from_dir() <moxel.html#moxel.utils.voxels_from_dir>`_

In all cases, :func:`Grid.calculate` is used under the hood to calculate the
energy voxels.

Functions :func:`voxels_from_file`, :func:`voxels_from_files`,
:func:`voxels_from_dir` are just wrappers. To better understand how to use
these functions check: :ref:`documentation`.

Visualization
^^^^^^^^^^^^^^^

.. code-block::

   >>> from moxel.utils import plot_voxels # Omitting .utils also works.
   >>> import matplotlib.pyplot as plt
   >>> fig = plot_voxels(voxels, cmap='coolwarm')
   >>> plt.show()

.. image:: images/plot_voxels.png
    :align: center
    :scale: 70%

Since ``voxels`` is just a ``np.array`` check also `plotly
<https://plotly.com/python/3d-volume-plots/>`_ and `pyvista
<https://docs.pyvista.org/version/stable/examples/02-plot/volume.html>`_.
