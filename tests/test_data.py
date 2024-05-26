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
from moxel.data import prepare_data
from moxel.utils import load_json


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
    ...


if __name__ == '__main__':
    unittest.main()
