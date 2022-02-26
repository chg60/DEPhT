import argparse
import json
import pathlib
import sys

from Bio import SeqIO, AlignIO

from depht.functions import multiprocess
from depht_utils.functions import clustalo

from depht_utils.data.defaults import HHSUITEDB_DEFAULTS as DEFAULTS


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_curate_functions(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)
    parser.add_argument("index_file", type=pathlib.Path)
    parser.add_argument("functions_config", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-np", "--cpus", type=int)

    parser.add_argument("-all", "--accept_all", action="store_true")
    parser.add_argument("-ms", "--min_size", type=int)

    parser.set_defaults(**DEFAULTS)
    args = parser.parse_args(unparsed_args)
    return args


def curate_functions(fasta_dir, index_file, functions_file, output_dir,
                     accept_all=False, verbose=False, cores=1,
                     min_size=10):
    cluster_functions = read_cluster_function_index_file(index_file)
    accept_list, ignore_list = read_functions_config_file(functions_file)

    function_to_cluster_map = map_function_to_cluster(cluster_functions)

    curated_clusters = []
    for cluster_function, clusters in function_to_cluster_map.items():
        ignore = False
        for ignored_function in ignore_list:
            if ignored_function.lower() in cluster_function.lower():
                ignore = True
                break

        if ignore:
            continue

        if not accept_all:
            for accepted_function in accept_list:
                if accepted_function.lower() in cluster_function.lower():
                    curated_clusters += clusters
                    break

            continue

        curated_clusters += clusters

    output_dir.mkdir(exist_ok=True, parents=True)
    work_items = []
    for cluster in curated_clusters:
        cluster_path = fasta_dir.joinpath(cluster)
        outpath = output_dir.joinpath(cluster)

        records = [record for record in SeqIO.parse(cluster_path, "fasta")]

        if len(records) < min_size:
            continue

        work_items.append((records, outpath, cluster_functions[cluster]))

    multiprocess.parallelize(work_items, cores, dump_named_alignment,
                             verbose=verbose)


def dump_named_alignment(records, fasta_path, name):
    with fasta_path.open(mode="w") as filehandle:
        SeqIO.write(records, filehandle, "fasta")

    clustalo.clustalo(fasta_path, fasta_path)

    aln_records = [record for record in AlignIO.read(fasta_path, "fasta")]

    with fasta_path.open(mode="w") as filehandle:
        filehandle.write("".join(["#", name, "\n"]))

        SeqIO.write(aln_records, filehandle, "fasta")


def map_function_to_cluster(cluster_functions):
    function_to_cluster_map = {}
    for cluster_name, cluster_function in cluster_functions.items():
        mapped_clusters = function_to_cluster_map.get(cluster_function, list())
        mapped_clusters.append(cluster_name)
        function_to_cluster_map[cluster_function] = mapped_clusters

    return function_to_cluster_map


def read_functions_config_file(functions_file):
    accepted_functions = []
    with functions_file.open(mode="r") as filehandle:
        functions_config = json.load(filehandle)

    accepted_functions = functions_config.get("LIKE", list())
    ignored_functions = functions_config.get("NOT LIKE", list())

    return accepted_functions, ignored_functions


def read_cluster_function_index_file(index_file):
    cluster_functions = {}
    with index_file.open(mode="r") as filehandle:
        lines = filehandle.readlines()
        for line in lines:
            split_line = line.split("\t")
            cluster_functions[split_line[0]] = split_line[1]

    return cluster_functions


def main(unparsed_args):
    args = parse_curate_functions(unparsed_args)
    curate_functions(args.input_dir, args.index_file, args.functions_config,
                     args.output_dir, cores=args.cpus, verbose=args.verbose,
                     accept_all=args.accept_all, min_size=args.min_size)


if __name__ == "__main__":
    main(sys.argv[1:])
