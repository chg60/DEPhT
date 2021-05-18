from pathlib import Path
from Bio import SeqIO
import argparse
import sys

from blastn import locate_subsequence

# The magic numbers
CUT_OFF = 0.5
CORRECT = 1
NOT_CORRECT = 0


def locate_seq(child_dir, parent_dir, csv_file):

    ends = {}

    for child_seq_path in child_dir.iterdir():

        child_base_name = child_seq_path.stem
        parent_filename = str(child_base_name) + ".fasta"
        parent_seq_path = Path(parent_dir/parent_filename)
        ends[child_base_name] = locate_subsequence(child_seq_path,
                                                   parent_seq_path,
                                                   csv_file)
    return ends


def get_lengths(manual_dir):

    lengths = {}

    for proph in dir.iterdir():
        with open(proph, "r") as fh:
            record = SeqIO.read(fh, "fasta")
            length = len(record)

            lengths[proph.stem] = length

    return lengths


def get_dir():

    """
    Parse argument list
    :return: Working directory where files are located
    :rtype: Path
    """

    parser = argparse.ArgumentParser(description=__doc__)
    dir_help = "Path to directory"

    parser.add_argument("dir", type=Path, help=dir_help)
    args = parser.parse_args()

    return args.dir


def compare_ends(ends_tuple, length):
    intersection = float(ends_tuple[1] - ends_tuple[0])
    # overlap = intersection/length  # Binary classifier

    # return overlap > CUT_OFF

    return (intersection/length) > CUT_OFF


def main():
    # Test variables
    '''
    p_tuple = (546, 1080)   # testing input
    r_tuple = (465, 1070)   # comparing against
    length_r = float(r_tuple[1] - r_tuple[0])
    p_set = set()
    r_set = set()

    # Making a set for the return tuple
    for i in range(p_tuple[0], p_tuple[1]+1):
        p_set.add(i)

    # Making a set for an existing tuple
    for i in range(r_tuple[0], r_tuple[1]+1):
        r_set.add(i)

    inter_set = p_set.intersection(r_set)

    # Overlapping score
    overlap = len(inter_set)/length_r

    # Binary classifier
    if overlap > CUT_OFF:
        return CORRECT
    return NOT_CORRECT

    dir_path = get_dir()
    results = {}

    for parent_dir in sorted(dir_path.iterdir()):
        lengths = get_lengths(parent_dir)
        results[parent_dir.name] = {}
        for child_dir in sorted(parent_dir.iterdir()):
            for proph in sorted(child_dir.iterdir()):
                csv_file = Path(child_dir_dir/"csv_file.csv")
                ends[child_dir.name] = {}
                ends[child_dir.name][proph] = locate_seq(child_dir,
                                              parent_dir, csv_file)

    for child_dir, proph in ends[c][p]:
        results[proph] = compare_ends(ends[proph], lengths[proph])

    print(results)
    '''

    main_dir = get_dir()

    if main_dir.is_dir() is False:
        print("Invalid directory")
        sys.exit(1)

    testing_data = []
    manual_data = []

    for parent_dir in sorted(main_dir.iterdir()):

        # for each directory in the parent directory - need to get results
        for child_dir in sorted(dir.iterdir()):
            if child_dir.is_dir() is False:     # parent sequence FASTA file
                continue

            csv_file = child_dir/"cvs_file.csv"

            ends = locate_seq(child_dir, parent_dir, csv_file)

            if child_dir.name == "Haley":
                manual_data.append(ends)
            else:
                testing_data.append(ends)

    print(manual_data)
    print(testing_data)


if __name__ == '__main__':
    main()
