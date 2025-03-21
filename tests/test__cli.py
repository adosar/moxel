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
import unittest
import tempfile


class TestCLI(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory(dir='/tmp')
        self.cif_dirname = 'tests/CIFs'
        self.out_pathname = os.path.join(self.tempdir.name, 'voxels_data')

    def run_command(self, command):
        exit_status = os.system(command)
        if exit_status != 0:
            raise RuntimeError(command)

    def test_cli_run(self):
        self.run_command(
                'moxel --config tests/config.yaml'
                f' {self.cif_dirname}'
                f' {self.out_pathname}'
                ' --grid_size 3'
                )

    def tearDown(self):
        self.tempdir.cleanup()


if __name__ == '__main__':
    unittest.main()
