|:chart_with_upwards_trend:| Changes
====================================

Version 1.0.0
-------------

.. attention::
    Required Python version changed from ``3.10`` to ``3.12``.


.. versionchanged:: 1.0.0

    * :func:`moxel.utils.voxels_from_files` and :func:`moxel.utils.voxels_from_dir`
      
        Now they store the names of the materials, so users don't need to do it.

.. versionremoved:: 1.0.0

    * :func:`moxel.utils.batch_clean_and_merge`

.. versionadded:: 1.0.0

    * :func:`moxel.utils.batch_clean`
    * :func:`moxel.visualize.plot_voxels_pv` for faster visualization.
    * Optional parameter ``n_jobs`` for specifying number of cores.
