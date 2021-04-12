from Bio import SeqIO
from pathlib import Path

import argparse
import sys
# import time


# Global variables
CUT_OFF = 0.5


# Will convert the file to a biopython object
def open_file(filename):

    """
    Converts the FASTA formatted file into a biopython object
    :param: filename
    :type: Path
    :return: SeqRecord object generated from the FASTA file
    :rtype: Bio.SeqRecord.SeqRecord
    """

    # print(filename)

    # Open file

    '''
    with open(filename, "r") as fh:
        # Convert file to a biopython seq object
        # use SeqIO.parse
        # only take in longest sequence!
        record = SeqIO.read(fh, "fasta")
    return record
    '''

    contigs = [record for record in SeqIO.parse(filename, "fasta")]
    contig_lengths = [len(record) for record in contigs]
    longest_contig = contigs[contig_lengths.index(max(contig_lengths))]

    return longest_contig


def find_sub(parent, other):

    """
    Locates the prophage region as a subsequence of the bacterial genome
    :param parent: Parent sequence
    :type parent: Bio.SeqRecord.SeqRecord
    :param other: Subsequence being located in the parent
    :type other: Bio.SeqRecord.SeqRecord
    :return: left end of the subsequence
    :rtype: int
    """

    # Find the beginning of the subsequence
    left_coor = parent.seq.find(other.seq)
    # Checks reverse complement
    if left_coor == -1:
        reverse = other.reverse_complement()
        left_coor = parent.seq.find(reverse.seq)
        # if nothing is returned an error exists
        if left_coor == -1:
            return None
    return left_coor


def compare_ends(manual_coor, test_coor):
    """
    Compares the manually annotated ends with the ends from testing software
    :param manual_coor: Coordinates from manually annotates genomes
    :param test_coor: Coordinates from testing software
    :return: True or False based on whether the overlap is over 50% or not
    :rtype: bool
    """
    manual_positions = set()
    test_positions = set()

    for i in range(manual_coor[0], manual_coor[1]+1):
        manual_positions.add(i)

    for i in range(test_coor[0], test_coor[1]+1):
        test_positions.add(i)

    inter_set = manual_positions.intersection(test_positions)

    # Overlapping score
    overlap = len(inter_set)/(manual_coor[1]-manual_coor[0])

    # return True/False
    return overlap > CUT_OFF


def get_dir():
    """
    Parses the argument list
    :return: Working directory where the sequences are located
    :rtype: Path
    """
    parser = argparse.ArgumentParser(description=__doc__)
    dir_help = "Path to working directory"

    parser.add_argument("dir", type=Path, help=dir_help)
    args = parser.parse_args()

    return args.dir


def get_parent_data(name, parent_dir):
    """
    Gets the ends of the parent sequence and the sequence contained in
    the parent FASTA
    :param name: Name of the parent genome
    :type name: str
    :param parent_dir: Directory where the parent FASTA is stored
    :return: List of two elements - tuple of ends and the parent sequence
    :rtype: list
    """
    parent_file = name + ".fasta"
    filepath = parent_dir/parent_file
    parent_seq = open_file(filepath)
    parent_tuple = (1, len(parent_seq))
    return [parent_tuple, parent_seq]


def get_child_ends(prophage_filepath, parent_seq):
    """
    Returns the ends of the prophage region as a tuple
    :param prophage_filepath: Filepath to the hypothetical prophage FASTA
    :type prophage_filepath: Path
    :param parent_seq: Sequence ofthe parent genome
    :return: Tuple of prophage ends
    :rtype: tuple
    """

    # change so that it reads all the prophages in that folder
    child_seq = open_file(prophage_filepath)
    start = find_sub(parent_seq, child_seq)

    if start is None:
        return None

    child_tuple = (start, start+len(child_seq))
    return child_tuple


def get_child_data(child_dir, parent_seq):

    """
    Returns a dictionary of the ends tuple with prophage name as the key
    :param child_dir: Child directory with hypothetical prophages
    :type child_dir: pathlib.Path
    :param parent_seq: Bacterial sequence
    :type parent_seq: SeqRecord object
    :return: Dictionary with prophage name keys and ends tuple values
    :rtype: dict
    """

    ends = {}

    for prophage in sorted(child_dir.iterdir()):
        prophage_name = prophage.stem
        ends[prophage_name] = get_child_ends(prophage, parent_seq)

    return ends


def print_data(manual_data, testing_data):
    """
    Prints the Testing Data
    :param manual_data: Data from the manually annotated genome
    :param testing_data: Data from testing software
    """

    print("\n----------------------------------------------------")
    print("Prophage\tEnds")
    print("----------------------------------------------------\n")
    for bac in manual_data:
        print(bac)
        print("Manually Annotated")
        for prophage in manual_data[bac]:
            print(f"{prophage}\t{manual_data[bac][prophage]}")

        for child_dir in testing_data[bac]:
            print("\n")
            print(child_dir)
            for prophage in testing_data[bac][child_dir]:
                print(f"{prophage}\t\t{testing_data[bac][child_dir][prophage]}"
                      )

        print("----------------------------------------------------\n")


def compile_data(main_dir):
    if main_dir.is_dir() is False:
        print("Invalid directory")
        sys.exit(1)

    testing_data = {}
    manual_data = {}

    for dir in sorted(main_dir.iterdir()):

        # data = {}
        # manual_coor = {}

        name = dir.name
        parent_data = get_parent_data(name, dir)
        # parent_ends = parent_data[0]
        parent_seq = parent_data[1]

        # dictionary with all the data for each phage
        testing_data[name] = {}
        manual_data[name] = {}

        # for each directory in the parent directory - need to get results
        for child_dir in sorted(dir.iterdir()):
            if child_dir.is_dir() is False:     # parent sequence FASTA file
                continue

            data_from = child_dir.name

            if data_from == "Haley":
                manual_coor = get_child_data(child_dir, parent_seq)
                manual_data[name] = manual_coor
            else:
                data = get_child_data(child_dir, parent_seq)
                testing_data[name][data_from] = data

    print_data(manual_data, testing_data)


def main():
    compile_data(get_dir())    # all the prophages are inside this directory


if __name__ == '__main__':
    main()
