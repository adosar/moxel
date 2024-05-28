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

import os
import json
import unittest
import tempfile
from itertools import combinations
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from moxel.utils import load_json
from moxel.data import prepare_data, VoxelDataset


class TestPrepareData(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory(dir='/tmp')
        self.filename = os.path.join(self.tempdir.name, 'clean_names.json')
        self.split_ratio = (6, 2, 2)
        self.all_indices = list(range(10))  # Dataset of size 10.

        # Create a dummy .json file.
        self.names = [f'name_{i}' for i in self.all_indices]
        with open(self.filename, 'w') as fhand:
            json.dump(self.names, fhand)

    def test_correct_split(self):
        # Split the data set.
        prepare_data(self.filename, split_ratio=self.split_ratio)

        # Load back the indices.
        train_idx, val_idx, test_idx = [
                load_json(os.path.join(self.tempdir.name, f'{mode}.json'))
                for mode in ['train', 'validation', 'test']
                ]

        # Their size must equal their split_ratio.
        for ds, size in zip([train_idx, val_idx, test_idx], self.split_ratio):
            self.assertEqual(len(ds), size)

        # Their pairwise intersections must be the empty set.
        for comb in combinations([train_idx, val_idx, test_idx], r=2):
            self.assertEqual(set(comb[0]) & set(comb[1]), set())

        # Their union must equal the pre-split indices.
        self.assertEqual(
                set(train_idx) | set(val_idx) | set(test_idx),
                set(self.all_indices)
                )

    def tearDown(self):
        self.tempdir.cleanup()


class TestVoxelsDataset(unittest.TestCase):
    def setUp(self):
        # The grid size is 5 and the train size is 4.
        self.path_to_indices = 'tests/toy_dataset/train.json'
        self.path_to_X = 'tests/toy_dataset/clean_voxels.npy'
        self.dummy_tfm_x = lambda x: x + 100
        self.ds_indices = load_json(self.path_to_indices)
        self.ds_voxels = np.load(self.path_to_X, mmap_mode='r')[:, None]
        self.batch_size = 2

        # For labeled dataset (optional).
        self.path_to_names = 'tests/toy_dataset/clean_names.json'
        self.path_to_Y = 'tests/toy_dataset/dummy.csv'
        self.index_col = 'id'
        self.labels = ['y_1', 'y_3']
        self.dummy_tfm_y = lambda x: x - 500

    def test_unlabeled_dataset(self):
        ds = VoxelDataset(
                path_to_indices=self.path_to_indices,
                path_to_X=self.path_to_X,
                transform_x=self.dummy_tfm_x,
                )

        dl = DataLoader(ds, batch_size=self.batch_size)

        # Check that it has correct size.
        self.assertEqual(len(ds), 4)

        # Check that it works with the dataloader.
        for x in dl:
            self.assertEqual(x.shape, (self.batch_size, 1, 5, 5, 5))  # Shape (B, C, D, H, W).
            self.assertIs(x.dtype, torch.float)

        # Check that transforms are correctly applied.
        for i in range(len(ds)):
            transformed_x = ds[i]
            raw_x = torch.tensor(self.ds_voxels[self.ds_indices[i]])
            self.assertTrue(torch.equal(transformed_x, self.dummy_tfm_x(raw_x)))

    def test_labeled_dataset(self):
        ds = VoxelDataset(
                path_to_indices=self.path_to_indices,
                path_to_X=self.path_to_X,
                path_to_names=self.path_to_names,
                path_to_Y=self.path_to_Y,
                index_col=self.index_col,
                labels=self.labels,
                transform_x=self.dummy_tfm_x,
                transform_y=self.dummy_tfm_y,
                )

        dl = DataLoader(ds, batch_size=self.batch_size)

        material_names = load_json(self.path_to_names)

        df = pd.read_csv(
                self.path_to_Y,
                index_col=self.index_col,
                usecols=[self.index_col, *self.labels],
                )

        # Check that it has correct size.
        self.assertEqual(len(ds), 4)

        # Check that it works with the dataloader.
        for x, y in dl:
            self.assertEqual(x.shape, (self.batch_size, 1, 5, 5, 5))  # Shape (B, C, D, H, W).
            self.assertIs(x.dtype, torch.float)
            self.assertEqual(y.shape, (self.batch_size, len(self.labels)))  # Shape (B, n_out).
            self.assertIs(y.dtype, torch.float)

        # Check that transforms are correctly applied.
        for i in range(len(ds)):
            idx = self.ds_indices[i]
            name = material_names[idx]

            transformed_x, transformed_y = ds[i]
            raw_x = torch.tensor(self.ds_voxels[idx])
            raw_y = torch.tensor(df.loc[name, self.labels].values)

            self.assertTrue(torch.equal(transformed_x, self.dummy_tfm_x(raw_x)))
            self.assertTrue(torch.equal(transformed_y, self.dummy_tfm_y(raw_y)))


if __name__ == '__main__':
    unittest.main()
