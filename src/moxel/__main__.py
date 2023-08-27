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
