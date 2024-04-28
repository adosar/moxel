|:pushpin:| Changes
===================

Version 1.0.0
-------------

.. attention::
    Required Python version changed from ``3.10`` to ``3.12``.


.. versionchanged:: 1.0.0

    * :func:`moxel.utils.voxels_from_files` and :func:`moxel.utils.voxels_from_dir`
      
        1. Now they store the names of the materials, so users don't need to do it.
        2. Parameter ``out_pathname`` now must be specified (no longer optional).

.. versionremoved:: 1.0.0

    * Easy imports, such as ``from moxel import Grid``.
    * :func:`moxel.utils.batch_clean_and_merge`

.. versionadded:: 1.0.0

    * :func:`moxel.utils.batch_clean`
    * :func:`moxel.visualize.plot_voxels_pv` for faster visualization.
    * Optional parameter ``n_jobs`` for specifying number of cores.
