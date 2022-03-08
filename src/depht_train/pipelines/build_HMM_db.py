"""Utility script to build an HHSuite3 HMM database of phams with
conserved phage functions.
"""
import argparse
import pathlib
import sys

from depht.data import GLOBAL_VARIABLES
from depht.functions.multiprocess import CPUS
from depht_train.functions import hhsuitedb

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
NAME = GLOBAL_VARIABLES["phage_homologs"]["essential_name"]


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_args(unparsed_args):
    """Parse commandline arguments.

    :param unparsed_args:
    :type unparsed_args: list[str]
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)

    parser.add_argument("-v", "--verbose",  action="store_true")
    parser.add_argument("-n", "--name", type=str)

    parser.add_argument("-np", "--cpus", type=int, default=CPUS)
    parser.add_argument("-mpi", "--use_mpi", action="store_true")

    args = parser.parse_args(unparsed_args)
    return args


def build_HMM_db(input_dir, output_dir, name=NAME,
                 cores=CPUS, use_mpi=False, verbose=False):
    output_dir.mkdir(exist_ok=True, parents=True)

    hhsuitedb.create_hhsuitedb(input_dir, output_dir, name, cores=cores,
                               use_mpi=use_mpi, verbose=verbose)


def main(unparsed_args):
    """Commandline entrypoint for this module."""
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])
    build_HMM_db(args.input_dir, args.output_dir, name=args.name,
                 cores=args.cpus, use_mpi=args.use_mpi,
                 verbose=args.verbose)


if __name__ == "__main__":
    main()
