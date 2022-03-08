"""Utility script to make a blastn database that DEPhT can use for
attB searches.
"""
import argparse
import pathlib
import sys

from Bio import SeqIO

from depht.data import GLOBAL_VARIABLES
from depht_train.functions import blastdb

NAME = GLOBAL_VARIABLES["reference_db"]["name"]


# MAIN FUNCTIONS
#  ----------------------------------------------------------------------------
def parse_args(unparsed_args):
    """Parse commandline arguments.

    :param unparsed_args:
    :type unparsed_args: list[str]
    """
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--name", type=str, default=NAME)

    args = parser.parse_args(unparsed_args)
    return args


def build_reference_db(input_dir, output_dir, name=NAME, verbose=False):
    """

    :param input_dir:
    :param output_dir:
    :param name:
    :param verbose:
    :return:
    """
    output_dir.mkdir(exist_ok=True, parents=True)

    cc_fasta_path = output_dir.joinpath(f"{name}.fasta")
    write_concatenated_fasta(input_dir, cc_fasta_path)

    blastdb.create_blastdb(cc_fasta_path, output_dir, name, verbose=verbose)


def write_concatenated_fasta(input_dir, cc_fasta_path):
    with cc_fasta_path.open(mode="w") as filehandle:
        for genome_file in input_dir.iterdir():
            if not genome_file.is_file():
                continue

            genome_records = SeqIO.parse(genome_file, "fasta")
            for genome_record in genome_records:
                filehandle.write(f">{genome_record.id}\n")
                filehandle.write(f"{str(genome_record.seq)}\n")


def main(unparsed_args=None):
    """Commandline entrypoint for this module."""
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])
    build_reference_db(args.input_dir, args.output_dir, name=args.name,
                       verbose=args.verbose)


if __name__ == "__main__":
    main()
