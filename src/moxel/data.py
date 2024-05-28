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
Write the docstring of the module.
"""

import os
import json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, random_split
from . utils import load_json


def prepare_data(source, split_ratio=(0.8, 0.1, 0.1), seed=1):
    r"""
    Split voxels into train, validation and test sets.

    .. warning::

        * You should use this function **after** :func:`utils.batch_clean`.
        * No directory is created by :func:`prepare_data`.
        * All ``.json`` files are stored under the directory containing ``source``.

    Each ``.json`` file stores the indices of ``clean_voxels.npy`` that will be
    used for training, validation and testing.

    Parameters
    ----------
    source: str
        Pathname to the file holding the names of the materials.
    split_ratio: sequence of size 3, default=(0.8, 0.1, 0.1)
        The sizes or fractions of splits to be produced.

        * ``train == split_ratio[0]``.
        * ``val == split_ratio[1]``.
        * ``test == split_ratio[2]``.

    seed : int, default=1
        Controls the randomness of the ``rng`` used for splitting.

    Examples
    --------

    Before the split::

        voxels_data
        ├──clean_voxels.npy
        └──clean_names.json

    >>> prepare_data('path/to/voxels_data/clean_names.json')  # doctest: SKIP

    After the split::

        voxels_data
        ├──clean_voxels.npy
        ├──clean_names.json
        ├──train.json
        ├──validation.json
        └──test.json
    """
    rng = torch.Generator().manual_seed(seed)
    path = Path(source).parent
    indices = range(len(load_json(source)))

    train, val, test = random_split(indices, split_ratio, generator=rng)

    for split, mode in zip((train, val, test), ('train', 'validation', 'test')):
        mode_indices = list(split)
        with open(os.path.join(path, f'{mode}.json'), 'w') as fhand:
            json.dump(mode_indices, fhand, indent=4)

    print('\033[32mData preparation completed!\033[0m')


class VoxelDataset(Dataset):
    r"""
    Dataset for voxels.

    .. _transforms: https://pytorch.org/tutorials/beginner/data_loading_tutorial.html#transforms

    .. tip::
        See `transforms`_ for implementing your own transforms.

    Parameters
    ----------
    path_to_indices: str
        Pathname to the ``.json`` file holding the indices of the voxels.
    path_to_X : str
        Pathname to the ``.npy`` file holding the voxels.
    path_to_names: str, optional
        Pathname to the ``.json`` file holding the names of the materials. No
        effect if ``path_to_Y=None``.
    path_to_Y : str, optional
        Pathname to the ``.csv`` file holding the labels of the voxels.
    index_col : str, optional
        Column name of the ``.csv`` file to be used as row labels. The names
        (values) under this column must follow the same naming scheme as in
        ``clean_names.json``.
    labels : list, optional
        List containing the column names of the ``.csv`` to be used as labels.
        No effect if ``path_to_Y=None``.
    transform_x : callable, optional
        Transforms applied to input.
    transform_y : callable, optional
        Transforms applied to output. No effect if ``path_to_Y=None``.
    """
    def __init__(
            self, path_to_indices, path_to_X,
            path_to_names=None, path_to_Y=None,
            index_col=None, labels=None,
            transform_x=None, transform_y=None,
            ):

        if (labels is not None) and (type(labels) is not list):
            raise ValueError('labels must be a list!')

        self.path_to_X = path_to_X
        self.path_to_Y = path_to_Y
        self.labels = labels
        self.index_col = index_col
        self.transform_x = transform_x
        self.transform_y = transform_y

        self._voxel_indices = load_json(path_to_indices)

        if path_to_Y is not None:
            self._voxel_names = load_json(path_to_names)

        self.X = None
        self.Y = None

    def __len__(self):
        return len(self._voxel_indices)

    def __getitem__(self, idx):
        if self.X is None:
            # Load and add a channel dimension to voxels array.
            self.X = np.load(self.path_to_X, mmap_mode='r')[:, None]

        voxel_idx = self._voxel_indices[idx]
        sample_x = torch.tensor(self.X[voxel_idx], dtype=torch.float)

        if self.transform_x is not None:
            sample_x = self.transform_x(sample_x)

        # Only for labeled datasets.
        if self.Y is None and self.path_to_Y is not None:
            self.Y = pd.read_csv(
                    self.path_to_Y,
                    index_col=self.index_col,
                    usecols=[*self.labels, self.index_col],
                    )

        if self.Y is not None:
            name = self._voxel_names[voxel_idx]
            sample_y = torch.tensor(self.Y.loc[name].values, dtype=torch.float)

            if self.transform_y is not None:
                sample_y = self.transform_y(sample_y)

            return sample_x, sample_y

        return sample_x
