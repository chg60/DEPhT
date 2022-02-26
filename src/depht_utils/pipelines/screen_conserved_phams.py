import argparse
import binascii
import math
import pathlib
import sys

import bitarray
from Bio import SeqIO

from depht_utils.data.defaults import SHELL_DB_DEFAULTS as DEFAULTS
from depht_utils.functions.fileio import (
    read_gene_index_file,
    read_cluster_index_file,
    write_gene_hex_value_file)

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
NAME = DEFAULTS["name"]

BINARY_TO_HEX_PLACES = 4

REP_THRESHOLD = 0.75


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_screen_conserved_phams(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)
    parser.add_argument("gene_index", type=pathlib.Path)
    parser.add_argument("cluster_index", type=pathlib.Path)

    parser.add_argument("-rt", "--representation_threshold", type=float,
                        default=REP_THRESHOLD)

    args = parser.parse_args(unparsed_args)
    return args


def screen_conserved_phams(input_dir, output_dir, gene_index, cluster_index,
                           rep_threshold=REP_THRESHOLD, name=NAME):
    gene_data = read_gene_index_file(gene_index)
    
    cluster_data = read_cluster_index_file(cluster_index)
    record_cluster_map = get_record_cluster_map(cluster_data)

    gene_rep_bitarrays = [None] * len(gene_data)
    for input_file in input_dir.iterdir():
        if input_file.suffix != ".fasta":
            continue
        
        gene_records = [record for record in SeqIO.parse(input_file, "fasta")]

        cluster_rep_bitarray = get_cluster_rep_bitarray(
                                            gene_records, gene_data,
                                            cluster_data, record_cluster_map,
                                            rep_threshold)
        
        for gene_record in gene_records:
            gene_rep_bitarrays[int(gene_record.id)] = cluster_rep_bitarray


    hex_length = math.ceil(len(cluster_data) / BINARY_TO_HEX_PLACES)
    gene_rep_hex_values = [binascii.b2a_hex(bitarray.tobytes())
                           for bitarray in gene_rep_bitarrays]

    gene_hex_value_file = output_dir.joinpath(name).with_suffix(".pbv")
    write_gene_hex_value_file(gene_hex_value_file, gene_rep_hex_values)


def get_cluster_rep_bitarray(gene_records, gene_data, cluster_data,
                             record_cluster_map, rep_threshold):
    cluster_representation_sets = [set() for _ in range(len(cluster_data))]

    for record in gene_records:
        gene_data_dict = gene_data.get(record.id, None)

        if gene_data_dict is None:
            continue

        parent = gene_data_dict["parent"]

        cluster_num = record_cluster_map.get(parent, None)

        if cluster_num is None:
            continue

        cluster_representation_sets[cluster_num].add(parent)

    bit_array_size = len(cluster_data)
    cluster_rep_bitarray = bitarray.bitarray([0] * bit_array_size)

    for i, cluster_rep_set in enumerate(cluster_representation_sets):
        representation = len(cluster_rep_set) / len(cluster_data[i])

        if representation < rep_threshold:
            continue

        cluster_rep_bitarray[i] = 1

    return cluster_rep_bitarray


def get_record_cluster_map(cluster_data):
    record_cluster_map = {}

    for index, record_names in enumerate(cluster_data):
        for record_name in record_names:
            record_cluster_map[record_name] = index

    return record_cluster_map


def main(unparsed_args):
    args = parse_screen_conserved_phams(unparsed_args)
    screen_conserved_phams(args.input_dir, args.output_dir, args.gene_index,
                           args.cluster_index,
                           rep_threshold=args.representation_threshold)


if __name__ == "__main__":
    main(sys.argv[1:])
