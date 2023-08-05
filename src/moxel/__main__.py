import os
import json
import argparse
from pathlib import Path
from . utils import voxels_from_dir


def _transaction_summary():
    col_size, _ = os.get_terminal_size()
    gap = col_size // 6
    num_cifs = len([i for i in os.listdir(args.directory) if i.endswith('.cif')])
    kelvin = 'K'
    angstrom = 'Å'

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
        f'{f"{args.cubic_box}/{args.centroid}/{args.length}":>{gap}}',
        sep=''
        )
    print('\nReading from directory:')
    print(f'  \033[1;31m{args.directory}\033[m')
    print('\n(Over)writing to file:')
    print(f'  \033[1;31m{args.o}\033[m')
    print('\nTransaction Summary')
    print(col_size*"=")
    print(f'Calculate voxels for {num_cifs} CIFs\n')


parser = argparse.ArgumentParser(
        prog='moxel',
        description='Generate voxelized potential energy surface \
                for a list of crystalline materials. CIF files are \
                    retrieved from <directory>.',
        epilog='For more information check: https://github.com/adosar/moxel.git'
        )

parser.add_argument('directory')
parser.add_argument(
        '-n', metavar='grid size',
        help='Number of grid points along each dimension.\
        \033[1;1mDefault=25\033[m',
        default=25, type=int
        )
parser.add_argument(
        '-c', metavar='cutoff',
        help='Cutoff radius for the Lennard-Jones potential.\
        \033[1;1mDefault=10\033[m',
        default=10, type=float
        )
parser.add_argument(
        '-e', metavar='epsilon',
        help='Epsilon (ε/K) value for the probe atom.\
        \033[1;1mDefault=50\033[m',
        default=50, type=float
        )
parser.add_argument(
        '-s', metavar='sigma',
        help='Sigma (σ/Å) value for the probe atom.\
        \033[1;1mDefault=2.5\033[m',
        default=2.5, type=float
        )
parser.add_argument(
        '--cubic-box',
        help='Set the simulation box to cubic.',
        action='store_true'
        )
parser.add_argument(
        '--centroid',
        help='The center of the cubic box in Å. Takes effect only if\
        flag \033[1;1m--cubic-box\033[m is set.\033[1;1m Default=0\033[m',
        default=0, type=float
        )
parser.add_argument(
        '--length',
        help='The size of the cubic box in Å. Takes effect only if\
        flag \033[1;1m--cubic-box\033[m is set.\033[1;1m Default=30\033[m',
        default=30, type=float
        )
parser.add_argument(
        '-o', metavar='output',
        help='Ouput file name to store the generated voxels.\
        \033[1;1mDefault=./voxels.npy\033[m',
        default='./voxels.npy'
        )

if __name__ == '__main__':

    #Loading LJ parameters
    with open(f'{Path(__file__).parents[0]}/lj_params.json', 'r') as fhand:
        lj_params = json.load(fhand)

    args = parser.parse_args()
    print(args)
    _transaction_summary()

    inp = input('\nIs this ok[y/N]: ')

    if inp == 'y':
        voxels_from_dir(
            args.directory, args.n,
            args.c, args.e, args.s,
            args.cubic_box, args.centroid,
            args.length, args.o
            )
    else:
        print('Operation aborted.')
