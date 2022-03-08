"""Utility script to retrieve sequences from Genbank using Entrez.
"""
import argparse
import pathlib
import sys

from Bio import SeqIO

from depht_train.functions import entrez

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULTS = {"input_type": "acc", "output_type": "fasta"}

INPUT_TYPES = ["acc", "tax"]
OUTPUT_TYPES = ["fasta", "gb"]


def parse_args(unparsed_args):
    """Parse commandline arguments.

    :param unparsed_args:
    :type unparsed_args: list[str]
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("identifiers_file", type=pathlib.Path)
    parser.add_argument("outdir", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument("-it", "--input_type", choices=INPUT_TYPES)
    parser.add_argument("-ot", "--output_type", choices=OUTPUT_TYPES)

    parser.set_defaults(**DEFAULTS)

    return parser.parse_args(unparsed_args)
    

def execute_pull_sequences(identifiers_file, output_dir, input_type="acc",
                           output_type="gb", verbose=False):
    """

    :param identifiers_file:
    :param output_dir:
    :param input_type:
    :param output_type:
    :param verbose:
    :return:
    """
    if not identifiers_file.is_file():
        print(f"File {identifiers_file} could not be found.")
        sys.exit(1)

    # FOR TESTING PLEASE REMOVE IF SEEN IN PRODUCTION
    entrez.set_entrez_credentials(email="chg60@pitt.edu")

    ids = entrez.parse_identifiers_file(identifiers_file)
    if input_type == "tax":

        accessions = set()
        for tax_id in ids:
            if verbose:
                print("Searching for sequence accessions "
                      f"with TaxonID {tax_id} ...")
            
                for acc_id in entrez.esearch_taxa(tax_id):
                    accessions.add(acc_id)
    elif input_type == "acc":
        accessions = ids
    else:
        raise NotImplementedError(f"Input type {input_type} is not supported.")

    output_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("Retrieving sequences with sequence accessions...")
    records = entrez.get_records(list(accessions))

    if verbose:
        print("Writing sequences to file...")
    for record in records:
        output_path = output_dir.joinpath(".".join([record.id, output_type]))

        with output_path.open(mode="w") as filehandle:
            SeqIO.write([record], filehandle, output_type)


def main(unparsed_args=None):
    """Commandline entrypoint to this module"""
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])

    execute_pull_sequences(args.identifiers_file, args.outdir,
                           input_type=args.input_type,
                           output_type=args.output_type,
                           verbose=args.verbose)


if __name__ == "__main__":
    main()
