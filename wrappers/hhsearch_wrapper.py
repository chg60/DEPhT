from pathlib import Path

from multiprocessing import cpu_count

from pdm_utils.functions.parallelize import *

import argparse
from subprocess import Popen, PIPE

def wrapper(input_file, output_dir, database, cutoff):

    # take in three arguments instead

    """
    Wrapper for HHSearch
    :param args: list of parsed arguments
    :type: list
    """

    output_file = Path(output_dir/(input_file.stem + ".hhr"))

    print(output_file)

    output = open(output_file, "w")

    # run hhsearch
    # subprocess.run(["hhsearch", "-i", input_file, "-d", database, "-e", cutoff])

    # return output file instead of out

    # with Popen(["hhsearch", "-i", input_file, "-d", database, "-e", cutoff], stdout=output, stderr=output) as p:
        # out = p.stdout.read().decode("utf-8")
        #err = p.stderr.read().decode("utf-8")

    Popen(["hhsearch", "-i", input_file, "-d", database, "-e", cutoff], stdout=output, stderr=output) 

    return output_file


def check_output_dir(output_dir):
    if output_dir.is_dir() is False:
        output_dir.mkdir()

    return output_dir


def create_job_queue(input_dir, output_dir, database, cutoff):

    jobs = []

    for filepath in input_dir.iterdir():
        if filepath.suffix == ".fasta":
            jobs.append((filepath, output_dir, database, cutoff))

    return jobs

def parse_args():

    """
    Parse arguments for the wrapper
    :returns: list containing parsed arguments
    """

    # ask for ouput, make sure it exists, stem of input for each - .hhr extension for output file
    # create if it doesn't
    # run subprocess with -o flag

    input_help = "Path to the input directory"
    output_help = "Path to the output directory"
    db_help = "Path to the database being searched against"
    cutoff_help = "E-value cutoff for results"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help=input_help)
    parser.add_argument("output_dir", type=Path, help=output_help)
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

    output_dir = check_output_dir(args.output_dir)

    jobs = create_job_queue(args.input_dir, output_dir, args.database, args.cutoff)

    results = parallelize(jobs, cpu_count()//2, wrapper) # list of output filepaths
    print(results)



if __name__ == "__main__":
    main()
