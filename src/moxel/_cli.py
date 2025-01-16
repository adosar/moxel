# This file is part of MOXελ.
# Copyright (C) 2023-2024 Antonios P. Sarikas

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

r"""
This module provides helper functions for the CLI.
"""

import fire
from . utils import voxels_from_dir


def moxel_fire():
    r"""
    CLI for the moxel package.
    """
    fire.Fire(voxels_from_dir)
