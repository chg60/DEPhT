"""Utility script to recut FASTA sequences to start at a new base"""

import sys
import argparse
import pathlib

from depht.functions.fasta import parse_fasta, write_fasta


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("infile", type=pathlib.Path, help="single-sequence FASTA file to recut")
    p.add_argument("start", type=int, help="which coordinate should become bp 1?")
    p.add_argument("outfile", type=pathlib.Path, help="output filepath")
    return p.parse_args()


def main():
    namespace = parse_args()
    infile = namespace.infile
    new_start = namespace.start
    outfile = namespace.outfile

    headers, sequences = parse_fasta(infile)
    if len(headers) != len(sequences):
        raise ValueError(f"{str(infile)} is not a valid FASTA file")
    if len(headers) > 1:
        raise ValueError(f"this script can only be run on single-entry FASTAs")

    sequence = sequences[0]
    sequence = sequence[new_start - 1:] + sequence[:new_start]
    sequences = [sequence]
    write_fasta(headers, sequences, outfile)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main()

