import argparse
import csv
import pathlib
import sys

from Bio import SeqIO

from prophicient_utils.data.defaults import HHSUITEDB_DEFAULTS
from prophicient_utils.functions.fileio import (
                                        read_cluster_table_file)


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
NAME = HHSUITEDB_DEFAULTS["name"]
TABLE = 11
DEFAULT_PRODUCT = HHSUITEDB_DEFAULTS["default_product"]


DEFAULTS = {"name": NAME, "table": TABLE, "default_product": DEFAULT_PRODUCT}


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_index_functions(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)
        
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--name", type=str)

    parser.add_argument("-cf", "--cluster_table", type=pathlib.Path,
                        default=None,)

    parser.set_defaults(**DEFAULTS)
    args = parser.parse_args(unparsed_args)
    return args


def index_functions(input_dir, output_dir,
                    name=DEFAULTS["name"], table=DEFAULTS["table"],
                    default_product=DEFAULTS["default_product"],
                    cluster_table=None):
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

                feature.qualifiers["note"] = [record.name]

                cds_features.append(feature)

    output_dir.mkdir(exist_ok=True, parents=True)

    fasta_file = output_dir.joinpath(".".join([name, "fasta"]))
    index_file = output_dir.joinpath(".".join([name, "pgi"]))

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

    if cluster_table is not None:
        clustered_records = get_clustered_records(cluster_table)

        cluster_file = output_dir.joinpath(".".join([name, "ci"]))

        with cluster_file.open(mode="w") as cluster_filehandle:
            for cluster_index, record_names in enumerate(clustered_records):
                cluster_filehandle.write("".join(
                                              [">", str(cluster_index), "\n"]))

                cluster_filehandle.write("".join(
                                              ["\0".join(record_names), "\n"]))


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


def main(unparsed_args):
    args = parse_index_functions(unparsed_args)
    index_functions(args.input_dir, args.output_dir, name=args.name,
                    cluster_table=args.cluster_table)


if __name__ == "__main__":
    main(sys.argv[1:])
