#!/usr/bin/python3

from pathlib import Path

import argparse
import os
import sys
import subprocess

'''
'''

current_working_directory = Path.cwd()



# def wrapper(input_path, output_path):
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

    formatted_args.extend(["-i", args.input_path, "-o", output_filepath])


    if args.training is not None:
        formatted_args.extend(["-t", args.training])
    if args.translation is not None:
        formatted_args.extend(["-a", args.translation])
    if args.nucleotide is not None:
        formatted_args.extend(["-d", args.nucleotide])
    if args.genes is not None:
        formatted_args.extend(["-s", args.genes])

    subprocess.run(formatted_args)


def main():

    """
    Parses the argument list from the command line and runs the wrapper
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
    parser.add_argument("-translation", type=Path, default=None, help=translation_help) # .txt
    parser.add_argument("-nucleotide", type=Path, default=None, help=nucleotide_help)
    parser.add_argument("-genes", type=Path, default=None, help=genes_help)

    args = parser.parse_args()

    # wrapper(args.input_path, args.output_path)
    wrapper(args)


if __name__ == "__main__":
    main()
