"""Pipeline to compile the complete set of databases and data structures
required to run DEPhT from bacterial and phage sequences provided.
"""
import argparse
import pathlib
import sys

from depht_utils import PACKAGE_DIR



def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    create_model()
    pass


def create_model():
    # Create a simple fasta-based database from the given bacterial sequences
    index_bacterial_sequences()

    # Phamerate bacterial sequences
    return
    # Create a simple fasta-based database from the given phage sequences
    index_phage_sequences() 


def index_phage_sequences():
    pass


def index_bacterial_sequences():
    pass


if __name__ == "__main__":
    unparsed_args = sys.argv
    if len(unparsed_args) <= 1:
        unparsed_args.append("-h")

    main(unparsed_args[1:])
