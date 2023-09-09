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
from . utils import voxels_from_dir
from . __init__ import __version__ as version


def _transaction_summary(args):
    col_size, _ = os.get_terminal_size()
    gap = col_size // 6
    num_cifs = len([i for i in os.listdir(args.directory) if i.endswith('.cif')])
    kelvin = 'K'
    angstrom = 'Å'
    out_info = args.o if args.o != None else './voxels.npy'

    print(col_size*"=")
    print(
        f'{"# CIFs":<{gap}}', f'{"Grid size":<{gap}}',
        f'{"Cutoff":<{gap}}', f'{"Epsilon":<{gap}}',
        f'{"Sigma":<{gap}}', f'{"Cubic-box":>{gap}}',
        sep=''
        )
    print(col_size*"=")
    print(
        f'{num_cifs:<{gap}}', f'{args.n:<{gap}}',
        f'{f"{args.c} {angstrom}":<{gap}}', f'{f"{args.e} {kelvin}":<{gap}}',
        f'{f"{args.s} {angstrom}":<{gap}}',
        f'{f"{args.cubic_box}/{args.length} Å":>{gap}}',
        sep=''
        )
    print('\nReading from directory:')
    print(f'  \033[1;31m{args.directory}\033[m')
    print('\n(Over)writing to file:')
    print(f'  \033[1;31m{out_info}\033[m')
    print('\nTransaction Summary')
    print(col_size*"=")
    print(f'Calculate voxels for {num_cifs} CIFs\n')


def _return_cli_parser():
    parser = argparse.ArgumentParser(
            prog='moxel-cli',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description='Generate energy voxels from a directory containig ``.cif`` files.',
            epilog='''A command line utility based on the MOXελ package.'''
            )

    parser.add_argument('--version', action='version', version=f'%(prog)s {version}')
    parser.add_argument('directory')
    parser.add_argument(
            '-n', metavar='grid_size',
            help='Number of grid points along each dimension.',
            default=25, type=int
            )
    parser.add_argument(
            '-c', metavar='cutoff',
            help='Cutoff radius (Å) for the LJ potential.',
            default=10, type=float
            )
    parser.add_argument(
            '-e', metavar='epsilon',
            help='Epsilon (ε/K) value for the probe atom.',
            default=50, type=float
            )
    parser.add_argument(
            '-s', metavar='sigma',
            help='Sigma (σ/Å) value for the probe atom.',
            default=2.5, type=float
            )
    parser.add_argument(
            '--cubic-box',
            help='Set the simulation box to cubic.',
            action='store_true'
            )
    parser.add_argument(
            '--length',
            help='The size of the cubic box in Å. Takes effect only if\
            flag **--cubic-box** is set.',
            default=30, type=float
            )
    parser.add_argument(
            '-o', metavar='output',
            help='Pathname to the file holding the voxels. If not specified,\
            voxels are stored in ``./voxels.npy``.',
            default=None
            )

    return parser
