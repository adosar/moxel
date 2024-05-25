r"""
Write the docstring of the module.
"""

import json
import torch
from torch.utils.data import Dataset, random_split
from . utils import load_json


def prepare_data(source, split_ratio=(0.8, 0.1, 0.1), seed=1):
    r"""
    Split a set of materials into train, validation and test sets.

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
    split_ratio: sequence, default=(0.8, 0.1, 0.1)
        The sizes or fractions of splits to be produced.
        * ``split_ratio[0] == train``.
        * ``split_ratio[1] == validation``.
        * ``split_ratio[2] == test``.
    seed : int, default=1
        Controls the randomness of the ``rng`` used for splitting.
    """
    rng = torch.Generator().manual_seed(seed)
    indices = range(len(load_json(source)))

    train, val, test = random_split(indices, split_ratio, generator=rng)

    for split, mode in zip((train, val, test), ('train', 'validation', 'test')):
        mode_indices = list(split)
        with open(os.path.join(path, f'{mode}.json'), 'w') as fhand:
            json.dump(mode_indices, fhand, indent=4)

    print('\033[32mData preparation completed!\033[0m')


class VoxelsDataset(Dataset):
    ...
