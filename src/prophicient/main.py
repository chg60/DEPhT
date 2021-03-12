"""
Prophicient scans bacterial genomes looking for prophages. Regions
identified as prophage candidates are further scrutinized, and
attachment sites identified as accurately as possible before
prophage extraction and generating the final report.
"""

import sys
import argparse
import pathlib

from Bio import SeqIO

from src.prophicient.functions.multiprocess import CPUS
from src.prophicient.functions.wrapper_basic import autoannotate
from src.prophicient.functions.prefilter import prefilter_genome


def parse_args(arguments):
    """
    Parse command line arguments
    :param arguments:  command line arguments that program was invoked with
    :type arguments: list
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", type=pathlib.Path,
                        help="FASTA file containing nucleotide sequence to scan for prophages")
    parser.add_argument("outdir", type=pathlib.Path,
                        help="path where output files can be written")
    parser.add_argument("--cpus", type=int, default=CPUS,
                        help=f"number of processors to use [default: {CPUS}]")
    return parser.parse_args(arguments)


def main(arguments):
    args = parse_args(arguments)

    # Verify that the input filepath is valid
    infile = args.infile
    if not infile.is_file():
        print(f"'{str(infile)}' is not a valid input file - exiting")
        sys.exit(1)

    # If the indicated output directory does not exist, make it
    outdir = args.outdir
    if not outdir.is_dir():
        print(f"'{str(outdir)}' does not exist - creating it")
        outdir.mkdir(parents=True)

    # Parse the input FASTA file
    records = list()
    with infile.open("r") as fasta_reader:
        record_iterator = SeqIO.parse(fasta_reader, "fasta")
        for record in record_iterator:
            records.append(record)

    # Create annotation outdir and auto-annotate
    annotate_dir = outdir.joinpath("prodigal")
    if not annotate_dir.is_dir():
        annotate_dir.mkdir(parents=False)
    prodigal_out = autoannotate(infile, annotate_dir)

    # TODO: stitch prodigal genes into records


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main(sys.argv[1:])
