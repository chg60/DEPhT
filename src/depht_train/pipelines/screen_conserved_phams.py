"""Utility script to identify conserved phams.
"""
import argparse
import binascii
import pathlib
import sys

import bitarray
from Bio import SeqIO

from depht.data import GLOBAL_VARIABLES
from depht_train.data import PARAMETERS
from depht_train.functions.fileio import (
    read_gene_index_file,
    read_cluster_index_file,
    write_gene_hex_value_file)

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
NAME = GLOBAL_VARIABLES["bacterial_sequences"]["name"]

REP_THRESHOLD = PARAMETERS["shell_db"]["rep_threshold"]

BINARY_TO_HEX_PLACES = 4


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_args(unparsed_args):
    """Parse commandline arguments.

    :param unparsed_args:
    :type unparsed_args: list[str]
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)
    parser.add_argument("gene_index", type=pathlib.Path)
    parser.add_argument("cluster_index", type=pathlib.Path)

    parser.add_argument("-n", "--name", type=str, default=NAME)
    parser.add_argument("-rt", "--representation_threshold", type=float,
                        default=REP_THRESHOLD)

    return parser.parse_args(unparsed_args)


def screen_conserved_phams(input_dir, output_dir, gene_index, cluster_index,
                           rep_threshold=REP_THRESHOLD, name=NAME):
    """Function to mark putative shell genome content simply by
    identifying subclade-specific conserved gene clusters.

    :param name:
    :param rep_threshold:
    :param input_dir: Path to dir containing protein cluster fasta-alignments
    :type input_dir: pathlib.Path
    :param output_dir:  Path to write conservation hexadecimal output
    :type output_dir: pathlib.Path
    :param gene_index: Path to gene index and relevant metadata file
    :type gene_index: pathlib.Path
    :param cluster_index: Path to a csv table mapping genomes to clusters.
    :type cluster_index: pathlib.Path
    """
    # Read in index and cluster data
    gene_data = read_gene_index_file(gene_index)
    cluster_data = read_cluster_index_file(cluster_index)
    record_cluster_map = get_record_cluster_map(cluster_data)

    # Initialize bitarray list data structure
    gene_rep_bitarrays = [None] * len(gene_data)
    for input_file in input_dir.iterdir():
        if input_file.suffix != ".fasta":
            continue

        gene_records = [record for record in SeqIO.parse(input_file, "fasta")]

        # Ascertain the conservation bitarray value from the members of
        # a protein sequence cluster
        cluster_rep_bitarray = get_cluster_rep_bitarray(
                                            gene_records, gene_data,
                                            cluster_data, record_cluster_map,
                                            rep_threshold)

        for gene_record in gene_records:
            gene_rep_bitarrays[int(gene_record.id)] = cluster_rep_bitarray

    # Convet the bitarrays to hexadecimal representations
    gene_rep_hex_values = [binascii.b2a_hex(bitarray.tobytes())
                           for bitarray in gene_rep_bitarrays]

    # Write the hexadecimal values in order found in the gene index
    gene_hex_value_file = output_dir.joinpath(name).with_suffix(".pbv")
    write_gene_hex_value_file(gene_hex_value_file, gene_rep_hex_values)


def get_cluster_rep_bitarray(gene_records, gene_data, cluster_data,
                             record_cluster_map, rep_threshold):
    """Function to calculate the bitarray value of a protein sequence cluster
    using the frequency of appearence within the clades of their respective
    hosts.

    :param gene_records: SeqRecord objects of the protein sequence cluster
    :type gene_records: list
    :param gene_data: Index and relevant metadata for each protein sequence
    :type gene_data: dict
    :param cluster_data: 2D-array of Sequence IDs grouped by clade membership
    :type cluster_data: list
    :param record_cluster_map: Sequence ID to cluster lookup map.
    :type record_cluster_map: dict
    :param rep_threshold: Cutoff for orthologous proteins to be shell content
    :type rep_threshold: float
    """
    cluster_representation_sets = [set() for _ in range(len(cluster_data))]

    # Create sets of all the host IDs for each protein sequence
    # in the sequence cluster
    for record in gene_records:
        gene_data_dict = gene_data.get(record.id, None)

        if gene_data_dict is None:
            continue

        parent = gene_data_dict["parent"]

        cluster_num = record_cluster_map.get(parent, None)

        if cluster_num is None:
            continue

        cluster_representation_sets[cluster_num].add(parent)

    # Intialize a bitarray to the size of the number of clades
    bit_array_size = len(cluster_data)
    cluster_rep_bitarray = bitarray.bitarray([0] * bit_array_size)

    # Flip bits depending on the frequency of appearence within
    # host sequence clades
    for i, cluster_rep_set in enumerate(cluster_representation_sets):
        representation = len(cluster_rep_set) / len(cluster_data[i])

        if representation < rep_threshold:
            continue

        cluster_rep_bitarray[i] = 1

    return cluster_rep_bitarray


def get_record_cluster_map(cluster_data):
    """

    :param cluster_data:
    :return:
    """
    record_cluster_map = {}

    for index, record_names in enumerate(cluster_data):
        for record_name in record_names:
            record_cluster_map[record_name] = index

    return record_cluster_map


def main(unparsed_args=None):
    """Commandline entrypoint to this module."""
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])
    screen_conserved_phams(args.input_dir, args.output_dir, args.gene_index,
                           args.cluster_index, name=args.name,
                           rep_threshold=args.representation_threshold)


if __name__ == "__main__":
    main()
