🚀 Tutorial
===========

As stated in :ref:`advantages`, all you need is a ``.cif`` file!

If you don't have one 👉 :download:`IRMOF-1.cif<down/IRMOF-1.cif>`.

Note that in the following examples, ``path/to/`` can be an absolute or relative
pathname.

Calculation for single material
-------------------------------

1. Functional interface:

    .. code-block::

        >>> from moxel.utils import voxels_from_file
        >>> voxels = voxels_from_file('path/to/IRMOF-1.cif', grid_size=25)

2. Object-oriented interface:

    .. code-block::

        >>> from moxel.utils import Grid
        >>> grid = Grid(grid_size=25)
        >>> grid.load_structure('path/to/IRMOF-1.cif')
        >>> grid.calculate()

.. code-block::

    >>> import numpy as np
    >>> np.array_equal(voxels, grid.voxels)  # A sanity check.
    True

.. tip::
    To visualize the energy voxels you can use PyVista::
    
        import pyvista as pv
        pl = pv.Plotter()
        pl.add_volume(voxels)
        pl.show()

    For interactive plots in Jupyter: ``pip install "pyvista[jupyter]"``

Of course, we are mostly interested in calculating voxels for multiple
materials. In this case, check:

* :func:`moxel.utils.voxels_from_files`
* :func:`moxel.utils.voxels_from_dir`

In all cases, :meth:`.Grid.calculate` is used under the hood to calculate the
voxels (all other functions are just wrappers). To better understand how to use
them refer to :ref:`documentation`.

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

Calculation for multiple materials
----------------------------------

Here, we examine how to generate energy voxels from a database, that can be
later used to train a ML algorithm (e.g. a CNN) or for manual feature extraction.

If you don't have a database 👉 :download:`CIFs.zip<down/CIFs.zip>`.

.. code-block:: console

   $ unzip path/to/CIFs.zip -d path/to/CIFs
   $ ls path/to/CIFs
   IRMOF-1.cif  ZnHBDC.cif  ZnMOF-74.cif

1. Create a directory to store voxels: 

    .. code-block:: console
        
        $ mkdir path/to/voxels_data

2. Calculate voxels and store them:

    .. tab-set::

        .. tab-item:: Python

            >>> from moxel.utils import voxels_from_dir
            >>> voxels_from_dir('path/to/CIFs/', 'path/to/voxels_data', grid_size=5)

        .. tab-item:: CLI

            .. code-block:: console

                $ moxel path/to/CIFs path/to/voxels_data/ --grid_size=5

            For reproducibility, it is recommended to use a configuration file::

                $ moxel --config=path/to/config.yaml

            You can generate and start customizing one as following::

                $ moxel --print_config > config.yaml

            You can also `override <https://jsonargparse.readthedocs.io/en/v4.37.0/#override-order>`_
            the values defined in the config::

                $ moxel --config=config.yaml --epsilon=100

            For more information: ``moxel --help``

        .. tab-item:: config.yaml

            .. literalinclude:: ../../tests/config.yaml
                :language: yaml


The voxels are stored as plain ``.npy`` files under ``path/to/voxels_data``:

    .. code-block:: console

        voxels_data/
        ├── IRMOF-1.npy
        ├── ZnHBDC.npy
        └── ZnMOF-74.npy
