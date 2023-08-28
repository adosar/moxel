.. highlight:: python

|:rocket:| Tutorial
===================

As stated in :ref:`advantages`, all you need is a ``.cif`` file!

If you don't have one |:point_right:| :download:`IRMOF-1.cif<down/IRMOF-1.cif>`.

Note that in the following examples, ``path/to/`` can be an absolute or relative pathname.

Calculation and visualization of voxels
---------------------------------------

Calculation
^^^^^^^^^^^

1. Functional interface:

    .. code-block::

        >>> from moxel.utils import voxels_from_file # Omitting .utils also works.
        >>> voxels = voxels_from_file('path/to/IRMOF-1.cif', grid_size=25)

2. Object-oriented interface:

    .. code-block::

        >>> from moxel.utils import Grid # Omitting .utils also works.
        >>> grid = Grid(grid_size=25)
        >>> grid.load_structure('path/to/IRMOF-1.cif')
        >>> grid.calculate()

.. code-block::

    >>> np.all(voxels == grid.voxels) # A sanity check.
    True

Of course, we are interested in calculating voxels from multiple files.
In this case, check:

* `voxels_from_files() <moxel.html#moxel.utils.voxels_from_files>`_
* `voxels_from_dir() <moxel.html#moxel.utils.voxels_from_dir>`_

In all cases, :func:`Grid.calculate` is used under the hood to calculate the
voxels. Functions :func:`voxels_from_file`, :func:`voxels_from_files`,
:func:`voxels_from_dir` are just wrappers. To better understand how to use these
functions check: :ref:`documentation`.

Visualization
^^^^^^^^^^^^^

.. code-block::

   >>> from moxel.utils import plot_voxels # Omitting .utils also works.
   >>> import matplotlib.pyplot as plt
   >>> import numpy as np
   >>> fill_pattern = np.tril(np.full(voxels.shape, True)) # Plot only the lower triangle.
   >>> fig = plot_voxels(voxels, fill_pattern=fill_pattern, cmap='coolwarm')
   >>> plt.show()

.. image:: images/plot_voxels.png
    :align: center
    :scale: 30%

Since ``voxels`` is just a ``np.array`` check also `plotly
<https://plotly.com/python/3d-volume-plots/>`_ and `pyvista
<https://docs.pyvista.org/version/stable/examples/02-plot/volume.html>`_.


Preparing voxels for a ML pipeline
----------------------------------

Here, we examine how to prepare clean ML inputs from a database, that can be
later used to train a ML algorithm (e.g. a CNN).

If you don't have a database |:point_right:| :download:`CIFs.zip<down/CIFs.zip>`.

.. code-block:: console

   $ unzip path/to/CIFs.zip -d path/to/CIFs
   $ ls path/to/CIFs
   corrupted_1.cif  corrupted_2.cif  IRMOF-1.cif  ZnHBDC.cif  ZnMOF-74.cif

Ideally, all ``.cif`` files should be processable. In this example, we cover the
general case where some ``.cif`` files (named as ``corrupted*``) can not be
processed.

1. Create a directory to store voxels: 

    .. code-block:: console
        
        $ mkdir path/to/batch 

2. Calculate voxels and store them:

   For this step you can also use the :ref:`cli`.

    .. tabs::

        .. code-tab:: python

            >>> from moxel.utils import voxels_from_dir # Omitting .utils also works.
            >>> voxels_from_dir('path/to/CIFs/', grid_size=10, out_name='path/to/batch/voxels.npy')

        .. code-tab:: console
            :caption: CLI

            $ python -m moxel -n 10 path/to/CIFs -o path/to/batch/voxels.npy

    Of course, it is necessary to know the indexing of the generated voxels.
    Under the hood, :func:`voxels_from_dir` uses
    ``sorted(os.listdir('path/to/CIFs/'))``, so we can just use a dictionary to
    keep track of the indexing::

        >>> import os, json
        >>> with open('path/to/batch/names.json', 'w') as fhand:
        ...     json.dump({'names': sorted(os.listdir('path/to/CIFs'))}, fhand, indent=4)


    .. warning::
        
        The directory structure::
            
            batch
            ├──voxels.npy
            └──names.json

        is necessary for :func:`batch_clean_and_merge`.

3. Clean the voxels:

    .. code-block::

        >>> from moxel.utils import batch_clean_and_merge # Omitting .utils also works.
        >>> exit_status = batch_clean_and_merge(['batch']) # You must pass a list!
        Missing voxels found! Cleaning...
        >>> exit_status
        1

    Lets check the contents of ``path/to/batch`` directory:

    .. code-block:: console
        
        $ ls path/to/batch
        clean_names.json  clean_voxels.npy  names.json  voxels.npy

    The file ``clean_names.json`` contains only the names of the processed
    materials:

    .. code-block:: console

        $ cat path/to/batch/clean_names.json
        {
            "names": [
                "IRMOF-1.cif",
                "ZnHBDC.cif",
                "ZnMOF-74.cif"
            ]
        }

    The file ``clean_voxels.npy`` contains only 3 samples:

    .. code-block::

        >>> import numpy as np
        >>> clean_voxels = np.load('path/to/batch/clean_voxels.npy', mmap_mode='r')
        >>> clean_voxels.shape
        (3, 10, 10, 10)

    
4. (optional) Remove ``voxels.npy`` and ``names.json``:

    .. code-block:: console

        $ rm path/to/batch/{voxels.npy,names.json}
