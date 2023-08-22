"""
MOXελ is a Python package for parallel calculation of energy voxels, with
emphasis on reticular chemistry.

.. note::
    Currently, interactions are modelled with Lennard-Jones potential.
"""

__author__ = 'Antonios P. Sarikas'
__version__ = '0.0.1'
__copyright__ = "Copyright (c) 2023 Antonios P. Sarikas"
__license__ = 'MIT'

from . utils import (
        Grid, voxels_from_file, voxels_from_files, voxels_from_dir,
        mic_scale_factors, batch_clean_and_merge, plot_voxels
        )
