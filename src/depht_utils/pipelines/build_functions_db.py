import argparse
import pathlib
import sys

from depht_utils.functions import hhsuitedb

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULTS = {"name": "Curated_Functions", "cpus": 1}


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_build_functions_db(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)

    parser.add_argument("-v", "--verbose",  action="store_true")
    parser.add_argument("-n", "--name", type=str)

    parser.add_argument("-np", "--cpus", type=int)
    parser.add_argument("-mpi", "--use_mpi", action="store_true")

    parser.set_defaults(**DEFAULTS)

    args = parser.parse_args(unparsed_args)
    return args


def build_functions_db(input_dir, output_dir, name=DEFAULTS["name"],
                       cores=1, use_mpi=False, verbose=False):
    output_dir.mkdir(exist_ok=True, parents=True)

    hhsuitedb.create_hhsuitedb(input_dir, output_dir, name, cores=cores,
                               use_mpi=use_mpi, verbose=verbose)


def main(unparsed_args):
    args = parse_build_functions_db(unparsed_args)
    build_functions_db(args.input_dir, args.output_dir, name=args.name,
                       cores=args.cpus, use_mpi=args.use_mpi,
                       verbose=args.verbose)


if __name__ == "__main__":
    main(sys.argv[1:])
