from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re

import argparse
import os
import sys
import subprocess

current_working_directory = Path.cwd()

PATTERN = re.compile(">\w+_(\d+) # (\d+) # (\d+) # (-?1) # ID=\w+;partial=(\d+);start_type=(\w+);rbs_motif=(.*|None);rbs_spacer=(.*bp|None);gc_cont=(\d.\d+)\s+")

def wrapper(args):

    """
    Wrapper function for prodigal
    :param args: list of arguments
    :type args: list
    """


    # the formatted arguments to run prodigal will be stored in a list
    formatted_args = ["prodigal"]

    # names the output file phage_name.output_format
    phage_name = args.input_path.stem

    output_format = ".gbk" #default

    output_filename = phage_name + output_format

    # creating the output filepath
    output_filepath = args.output_path/output_filename

    formatted_args.extend(["-i", args.input_path, "-a", output_filepath])


    if args.training is not None:
        formatted_args.extend(["-t", args.training])
    if args.nucleotide is not None:
        formatted_args.extend(["-d", args.nucleotide])
    if args.genes is not None:
        formatted_args.extend(["-s", args.genes])

    subprocess.run(formatted_args)

    return output_filepath


def parse_output_file(input_file, output_file):

    """
    Parses the output file from Prodigal
    :param output_file: output file from Prodigal
    :type output_file: Path
    """

    with open(input_file, "r") as fh:
        record = SeqIO.read(fh, "fasta")

    # Read the Prodigal output and parse it using PATTERN
    with open(output_file, "r") as fh:
        contents = "".join(fh.readlines())

    prodigal_genes = PATTERN.findall(contents)

    # Iterate over prodigal_genes (regex "hits" from Prodigal file...)
    for gene in prodigal_genes:

        gene_num = int(gene[0])
        start = int(gene[1])
        end = int(gene[2])
        strand = int(gene[3])
        partial = int(gene[4])
        start_codon = gene[5]
        rbs_type = gene[6]
        rbs_spacer = gene[7]
        gc_pct = gene[8]

        # Create SeqFeature from these data, and add it to record.features
        qualify={"note":{"partial": partial, "start_codon": start_codon,
                "rbs_type": rbs_type, "rbs_spacer": rbs_spacer,
                "gc content": gc_pct} }
        feature = SeqFeature(FeatureLocation(start, end),strand, qualifiers=qualify)
        record.features.append(feature)

def parse_args():

    """
    Parses the argument list from the command line and runs the wrapper
    :return args: list of parsed arguments
    """

    input_help = "Path to the input file"
    output_help = "Path to the output directory"
    training_help = "Use specified training file or create one if not present"
    translation_help = "Store the product of every gene to the specified file"
    nucleotide_help = "Store the nucleotide sequences in the specified file"
    genes_help = "Write all potential genes with scores to the selected file"


    # add help msgs
    parser = argparse.ArgumentParser(description=__doc__) # description = __doc__
    parser.add_argument("input_path", type=Path, help=input_help)
    parser.add_argument("output_path", type=Path, help=output_help)
    parser.add_argument("-training", type=Path, default=None, help=training_help)
    parser.add_argument("-nucleotide", type=Path, default=None, help=nucleotide_help)
    parser.add_argument("-genes", type=Path, default=None, help=genes_help)

    args = parser.parse_args()

    return args



def main():

    # args: input path, output path, modes
    args = parse_args()

    input_file = args.input_path

    # run the Prodigal wrapper
    output_file = wrapper(args)
    # parse the output file (.gbk)
    parse_output_file(input_file, output_file)


if __name__ == "__main__":
    main()
