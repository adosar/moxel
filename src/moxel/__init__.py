# This file is part of MOXελ.
# Copyright (C) 2023 Antonios P. Sarikas

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

"""
MOXελ is a Python package for parallel calculation of energy voxels, with
emphasis on reticular chemistry.

.. note::
    Currently, interactions are modelled with the Lennard-Jones (LJ) potential.
"""

__author__ = 'Antonios P. Sarikas'
__version__ = '0.0.0'
__copyright__ = "Copyright (c) 2023 Antonios P. Sarikas"
__license__ = 'GPL-3.0-only'

from . utils import (
        Grid, voxels_from_file, voxels_from_files, voxels_from_dir,
        mic_scale_factors, batch_clean_and_merge, plot_voxels
        )
