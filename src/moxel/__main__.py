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
import argparse
from pathlib import Path
from . cli import _transaction_summary, _return_cli_parser
from . utils import voxels_from_dir


if __name__ == '__main__':

    args = _return_cli_parser().parse_args()
    _transaction_summary(args)

    inp = input('\nIs this ok[y/N]: ')
    print('\n')

    if inp == 'y':
        voxels_from_dir(
            args.directory, args.n,
            args.c, args.e, args.s,
            args.cubic_box, args.length,
            args.o
            )
    else:
        print('Operation aborted.')
