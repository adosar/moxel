|:pushpin:| Changes
===================

Version 0.1.1
-------------
    
* Add performance enhancements proposed by :user:`fxcoudert`. See PR :pr:`2, 3`.
* Fix bug in :mod:`moxel.visualize` where plots didn't render.

Version 0.1.0
-------------

.. versionchanged:: 0.1.0

    * :func:`moxel.utils.voxels_from_files` and :func:`moxel.utils.voxels_from_dir`
      
        1. Now they store the names of the materials as a :class:`list`,
           so users don't need to do it.
        2. Parameter ``out_pathname`` now must be specified (no longer optional).

    * The usage of the CLI is now ``moxel <command>`` instead of ``python -m
      moxel``.

.. versionremoved:: 0.1.0

    * Easy imports, such as ``from moxel import Grid``.
    * :func:`moxel.utils.batch_clean_and_merge`

.. versionadded:: 0.1.0

    * :func:`moxel.utils.batch_clean`
    * :func:`moxel.visualize.plot_voxels_pv` for faster visualization.
    * Optional parameter ``n_jobs`` for specifying number of cores.
