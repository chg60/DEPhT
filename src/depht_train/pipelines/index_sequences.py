"""Utility script to index protein-coding gene sequences.
"""
import argparse
import pathlib
import sys

from Bio import SeqIO

from depht.data import GLOBAL_VARIABLES
from depht_train.functions.fileio import (
    read_cluster_table_file, write_cluster_file)

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
TABLE = 11
NAME = GLOBAL_VARIABLES["sequences"]["name"]
DEFAULT_PRODUCT = GLOBAL_VARIABLES["sequences"]["default_product"]


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_args(unparsed_args):
    """Parse commandline arguments

    :param unparsed_args: command line args
    :type unparsed_args: list[str]
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--name", type=str, default=NAME)

    parser.add_argument("-tt", "--translation_table", type=int,
                        default=TABLE)
    parser.add_argument("-cf", "--cluster_table", type=pathlib.Path,
                        default=None)

    args = parser.parse_args(unparsed_args)
    return args


def index_sequences(input_dir, output_dir,
                    name=NAME, table=TABLE,
                    default_product=DEFAULT_PRODUCT,
                    cluster_table=None):
    """Function to index protein coding features and related metadata.

    :param input_dir: Path to annotated GenBank-formatted genome sequences.
    :type input_dir: pathlib.Path
    :param output_dir: Path to write index output.
    :type output_dir: pathlib.Path
    :param name: Stem name of the outputted index file.
    :type name: str
    :param table: Translation table key to translate sequences, if necessary
    :type table: int
    :param default_product: Product to annotated unannotated sequences with.
    :type default_product: str
    :param cluster_table: Path to a csv table mapping genomes to clusters.
    :type cluster_table: pathlib.Path
    """
    cds_features = []
    for input_file in input_dir.iterdir():
        records = [record for record in SeqIO.parse(input_file, "gb")]

        for record in records:
            for index, feature in enumerate(record.features):
                if feature.type != "CDS":
                    continue

                # Checks to see if the CDS feature has an annotated translation
                translation = feature.qualifiers.get("translation")
                # If the feature lacks a translation qualifier, create one
                if not translation:
                    sequence = feature.extract(record.seq)
                    translation = sequence.translate(table=table, to_stop=True)

                    feature.qualifiers["translation"] = [translation]

                # Checks to see if the CDS feature has an annotated function
                product = feature.qualifiers.get("product")
                # If the feature lacks a function qualifier, create one
                if not product:
                    feature.qualifiers["product"] = [default_product]
                else:
                    # If the feature has a blank function qualifier, create one
                    if not product[0]:
                        feature.qualifiers["product"] = [default_product]

                # Checks to see if the CDS feature has an annotated locus tag
                locus_tag = feature.qualifiers.get("locus_tag")
                if not locus_tag:
                    feature.qualifiers["locus_tag"] = ["_".join(
                                                    [record.id, str(index+1)])]

                feature.qualifiers["note"] = [input_file.stem]

                cds_features.append(feature)

    output_dir.mkdir(exist_ok=True, parents=True)

    fasta_file = output_dir.joinpath(".".join([name, "fasta"]))
    index_file = output_dir.joinpath(".".join([name, "pgi"]))
    write_index_files(cds_features, index_file, fasta_file)

    cluster_file = None
    if cluster_table is not None:
        clustered_ids = get_clustered_records(cluster_table)

        cluster_file = output_dir.joinpath(".".join([name, "ci"]))
        write_cluster_file(clustered_ids, cluster_file)

    return fasta_file, index_file, cluster_file


def get_clustered_records(cluster_table_path):
    record_cluster_map = read_cluster_table_file(cluster_table_path)

    cluster_records_map = dict()
    for name, cluster in record_cluster_map.items():
        if cluster == "":
            cluster = None

        names = cluster_records_map.get(cluster, list())
        names.append(name)
        cluster_records_map[cluster] = names

    clustered_records = [names for c, names in cluster_records_map.items()
                         if c is not None]

    return clustered_records


def write_index_files(cds_features, index_file, fasta_file):
    """Function to write indexed sequences to file.

    :param cds_features: SeqFeature objects to be indexed.
    :type cds_features: list
    :param fasta_file: Path to write indexed sequences in fasta-format.
    :type fasta_file: pathlib.Path
    :param index_file: Path to write index and relevant metadata.
    :type index_file: pathlib.Path
    """
    with fasta_file.open(mode="w") as fasta_filehandle:
        with index_file.open(mode="w") as index_filehandle:
            for index, feature in enumerate(cds_features):
                locus_tag = feature.qualifiers["locus_tag"][0]
                translation = feature.qualifiers["translation"][0]
                product = feature.qualifiers["product"][0]
                parent = feature.qualifiers["note"][0]

                fasta_filehandle.write("".join([">", str(index), "\n"]))
                fasta_filehandle.write("".join([translation, "\n"]))

                index_filehandle.write("\t".join(
                                       [str(index), locus_tag, product,
                                        parent, "\n"]))


def main(unparsed_args=None):
    """Commandline entrypoint to this module."""
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])
    index_sequences(args.input_dir, args.output_dir, name=args.name,
                    table=args.translation_table,
                    cluster_table=args.cluster_table)


if __name__ == "__main__":
    main()
