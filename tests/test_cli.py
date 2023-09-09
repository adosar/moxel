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
Run tests from the project's root.

python -m unittest tests.<test_module>
"""

import os
import unittest
import tempfile


class TestMoxelCLI(unittest.TestCase):

    def test_help(self):
        exit_status = os.system('python -m src.moxel -h')
        self.assertEqual(exit_status, 0)

    def test_run(self):
        with tempfile.TemporaryDirectory() as dir_path:
            file_name = 'voxels.npy'
            exit_status = os.system(
                    f'python -m src.moxel tests/CIFs/foo\
                            -n 5 -o {dir_path}/{file_name}'
                    )
            self.assertEqual(exit_status, 0)


if __name__ == '__main__':
    unittest.main()
