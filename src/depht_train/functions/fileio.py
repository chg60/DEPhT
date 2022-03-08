import csv
import json


# READ FUNCTIONS
# -----------------------------------------------------------------------------

def read_gene_index_file(index_file):
    gene_index = {}
    with index_file.open(mode="r") as filehandle:
        lines = filehandle.readlines()
        for line in lines:
            split_line = line.split("\t")
            gene_index[split_line[0]] = {"locus_tag": split_line[1],
                                         "product": split_line[2],
                                         "parent": split_line[3]}

    return gene_index


def read_cluster_table_file(cluster_table_path):
    """Function to read in a cluster table csv

    :param cluster_table_path: Path to a csv table mapping genomes to clusters.
    :type cluster_table: pathlib.Path
    """
    # Use DictReader to read in all rows in the csv as dictionaries,
    # mapping headers to row values
    data_dicts = []
    with cluster_table_path.open(mode="r") as filehandle:
        filereader = csv.DictReader(filehandle)
        for data_dict in filereader:
            data_dicts.append(data_dict)

    # Create a lookup dictionary of genome names to their respective clusters
    record_cluster_map = dict()
    for data_dict in data_dicts:
        record_cluster_map[data_dict["Name"]] = data_dict["Cluster"]

    return record_cluster_map


def read_cluster_index_file(cluster_index):
    with cluster_index.open(mode="r") as filehandle:
        lines = filehandle.readlines()

    cluster_data = []

    line_index = 0
    while True:
        if line_index >= len(lines):
            break

        line = lines[line_index]
        if line.startswith(">"):
            line_index += 1

            if line_index >= len(lines):
                break

            cluster_members = [line.rstrip()
                               for line in lines[line_index].split("\0")]

            cluster_data.append(cluster_members)

        line_index += 1

    return cluster_data


def read_functions_config_file(functions_file):
    accepted_functions = []
    with functions_file.open(mode="r") as filehandle:
        functions_config = json.load(filehandle)

    accepted_functions = functions_config.get("LIKE", list())
    ignored_functions = functions_config.get("NOT LIKE", list())

    return accepted_functions, ignored_functions


# WRITE FUNCTIONS
# -----------------------------------------------------------------------------
def write_cluster_file(clustered_ids, cluster_file):
    with cluster_file.open(mode="w") as cluster_filehandle:
        for cluster_index, record_names in enumerate(clustered_ids):
            cluster_filehandle.write("".join(
                                          [">", str(cluster_index), "\n"]))

            cluster_filehandle.write("".join(
                                          ["\0".join(record_names), "\n"]))


def write_cluster_function_index_file(gene_cluster_data, index_file):
    with index_file.open(mode="w") as filehandle:
        for gene_cluster in gene_cluster_data:
            filehandle.write("\t".join(gene_cluster + ["\n"]))


def write_gene_hex_value_file(gene_hex_value_file, gene_rep_hex_values):
    with gene_hex_value_file.open(mode="wb") as filehandle:
        for hex_value in gene_rep_hex_values:
            filehandle.write(hex_value)
            filehandle.write(b"_")
