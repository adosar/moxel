r"""
Write the docstring of the module.
"""

import os
import json
from pathlib import Path
import torch
from torch.utils.data import Dataset, random_split
from . utils import load_json


def prepare_data(source, split_ratio=(0.8, 0.1, 0.1), seed=1):
    r"""
    Split a list of materials' names into train, validation and test sets.

    .. warning::
        * You should use this function **after** :func:`utils.batch_clean`.
        * No directory is created by :func:`prepare_data`. **All ``.json``
        files are stored under the directory containing ``source``**.

    Before the split::

        voxels_data
        ├──clean_voxels.npy
        └──clean_names.json

    After the split::

        voxels_data
        ├──clean_voxels.npy
        ├──clean_names.json
        ├──train.json
        ├──validation.json
        └──test.json

    Each ``.json`` file stores the indices of ``clean_voxels.npy`` that will be
    used for training, validation and testing.

    Parameters
    ----------
    source: str
        Pathname to the file holding the names of the materials
        (``clean_names.json``).
    split_ratio: array-like of shape (3,), default=(0.8, 0.1, 0.1)
        The sizes or fractions of splits to be produced.
        * ``split_ratio[0] == train``.
        * ``split_ratio[1] == validation``.
        * ``split_ratio[2] == test``.
    seed : int, default=1
        Controls the randomness of the ``rng`` used for splitting.
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

    Parameters
    ----------
    path_to_indices: str
        Pathname to the ``.json`` file holding the indices of the voxels.
    path_to_X : str
        Pathname to the ``.npy`` file holding the voxels.
    path_to_names: str, optional
        Pathname to the ``.json`` file holding the names of the materials. No
        effect if ``path_to_Y == None``.
    path_to_Y : str, optional
        Pathname to the ``.csv`` file holding the labels of the voxels.
    index_col : str, optional
        Column name of the ``.csv`` file to be used as row labels. The values
        (names) under this column must follow the same naming scheme as in
        ``clean_names.json``.
    labels : list, optional
        List containing the names of the properties to be predicted. No effect
        if ``path_to_Y == None``.
    transform_x : callable, optional
        Transforms applied to input. See `transforms`_ for implementing your own
        transforms.
    transform_y : callable, optional
        Transforms applied to output.  See `transforms`_ for implementing your
        own transforms. No effect if ``pcd_Y == None``.

        .. note::
            For example, if you want to perform classification, here you can
            pass the one-hot encoder (if the dataset is not already preprocessed).

    .. _transforms: https://pytorch.org/tutorials/beginner/data_loading_tutorial.html#transforms
    """
    def __init__(
            self, path_to_indices, path_to_X,
            path_to_names=None, path_to_Y=None,
            index_col=None, labels=None,
            transform_x=None, transform_y=None,
            ):

        if (labels is not None) and (type(labels) != list):
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

    @property
    def voxel_indices(self):
        return self._voxel_indices

    def __len__(self):
        return len(self.voxel_indices)

    def __getitem__(self, idx):
        # Account for np.load and multiprocessing.
        if self.X is None:
            self.X = np.load(self.path_to_X, mmap_mode='r')
        if self.Y is None and self.path_to_Y is not None:
            self.Y = pd.read_csv(
                    self.path_to_Y,
                    index_col=self.index_col,
                    usecols=[*self.labels, self.index_col],
                    )

        sample_x = self.X[name]

        if self.transform_x is not None:
            sample_x = self.transform_x(sample_x)

        # Only for labeled datasets.
        if self.Y is not None:
            name = self._voxel_names[idx]
            sample_y = self.Y.loc[name].values

            if self.transform_y is not None:
                sample_y = self.transform_y(sample_y)

            return (
                    torch.tensor(sample_x, dtype=torch.float),
                    torch.tensor(sample_y, dtype=torch.float)
                    )

        return torch.tensor(sample_x, dtype=torch.float)
