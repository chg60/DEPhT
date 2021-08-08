import argparse
import pathlib
import sys

from Bio import SeqIO

from depht_utils.functions import blastdb

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULTS = {"name": "Mycobacteria"}


# MAIN FUNCTIONS
#  ----------------------------------------------------------------------------
def parse_build_reference_db(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--name", type=str)
   
    parser.set_defaults(**DEFAULTS)

    args = parser.parse_args(unparsed_args)
    return args


def build_reference_db(input_dir, output_dir, name=DEFAULTS["name"],
                       verbose=False):
    output_dir.mkdir(exist_ok=True, parents=True)

    cc_fasta_path = output_dir.joinpath("references.fasta") 
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


def main(unparsed_args):
    args = parse_build_reference_db(unparsed_args)
    build_reference_db(args.input_dir, args.output_dir, name=args.name,
                       verbose=args.verbose)


if  __name__ == "__main__":
    main(sys.argv[1:])
