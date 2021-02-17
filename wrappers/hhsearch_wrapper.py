from pathlib import Path

import argparse
import subprocess

def wrapper(args):

    """
    Wrapper for HHSearch
    :param args: list of parsed arguments
    :type: list
    """


    input_file = args.input_file
    database = args.database
    cutoff = args.cutoff


    # run hhsearch
    subprocess.run(["hhsearch", "-i", input_file, "-d", database, "-e", cutoff])

def parse_args():

    """
    Parse arguments for the wrapper
    :returns: list containing parsed arguments
    """

    input_help = "Path to the input file"
    db_help = "Path to the database being searched against"
    cutoff_help = "E-value cutoff for results"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", type=Path, help=input_help)
    parser.add_argument("database", type=Path, help=db_help)
    parser.add_argument("-c", "--cutoff", type=str, default="0.001", help=cutoff_help)
    args = parser.parse_args()

    return args

def main():

    """
    Runs the wrapper
    """

    args = parse_args()

    print(args)
    wrapper(args)



if __name__ == "__main__":
    main()
