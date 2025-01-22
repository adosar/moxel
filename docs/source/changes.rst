📌 Changelog
============

Version 0.2.0
-------------

.. versionchanged:: 0.2.0

   * Storing scheme for energy voxels. They are now stored as (individual) plain
     ``.npy`` files (:issue:`8`).
    
.. versionremoved:: 0.2.0

   * :mod:`moxel.visualize` since PyVista already provides a simple way to
     visualize voxels (:issue:`9`).
   * :func:`~moxel.utils.load_json` which was needed by :func:`~moxel.utils.batch_clean`.
   * :func:`~moxel.utils.batch_clean` since it is no longer necessary to "clean"
     voxels due to new storing scheme.

Version 0.1.2
-------------

.. versionadded:: 0.1.2

    * :func:`moxel.utils.load_json`
    * Documentation for the CLI.

.. versionremoved:: 0.1.2

    * :func:`moxel.utils.get_names`

Version 0.1.1
-------------
    
* Performance enhancements proposed by :user:`fxcoudert`. See PR :pr:`2, 3`.
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

.. versionadded:: 0.1.0

    * :func:`moxel.utils.batch_clean`
    * :func:`moxel.visualize.plot_voxels_pv` for faster visualization.
    * Optional parameter ``n_jobs`` for specifying number of cores.

.. versionremoved:: 0.1.0

    * Easy imports, such as ``from moxel import Grid``.
    * :func:`moxel.utils.batch_clean_and_merge`

