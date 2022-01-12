import argparse
import pathlib
import sys

from Bio import SeqIO, AlignIO

from depht.functions import multiprocess
from depht_utils.functions import clustalo, fileio
from depht_utils.data.defaults import HHSUITEDB_DEFAULTS

# GLOBAL VARAIABLES
DEFAULTS = {"cpus": 1, "min_size": 10, "cutoff": 0.3}



# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_curate_functions(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)
    parser.add_argument("index_file", type=pathlib.Path)
    parser.add_argument("functions_config", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--name", type=str)
    parser.add_argument("-np", "--cpus", type=int)

    parser.add_argument("-c", "--cutoff", type=float)
    parser.add_argument("-all", "--accept_all", action="store_true")
    parser.add_argument("-ms", "--min_size", type=int)

    parser.set_defaults(**DEFAULTS)
    args = parser.parse_args(unparsed_args)
    return args


def curate_gene_clusters(fasta_dir, index_file, functions_file, output_dir,
                         accept_all=False, verbose=False, cores=1,
                         min_size=10, cutoff=DEFAULTS["cutoff"]):
    """Annotate and curate gene clusters, discarding those gene clusters
    with undesired or unclear consensus functions

    :param fasta_dir: Path to dir containing protein cluster fasta-alignments
    :type fasta_dir: pathlib.Path
    :param index_file: Path to gene index and relevant metadata file
    :type index_file: pathlib.Path
    :param functions_file: JSON parameters file for desired/undesired products
    :type functions_file:pathlib.Path
    :param output_dir: Path to directory to dump named/aligned sequence files
    :type output_dir: pathlib.Path
    :param cores: Number of processes to utilize
    :type cores: int
    :param verbose: Toggle pipeline print statements
    :type verbose: bool
    :param min_size: Minimum size threshold for protein clusters
    :type min_size: int
    :param cutoff: Frequency fraction rerquired for the consensus product.
    :type cutoff: float
    """
    # Annotate gene clusters based on consensus product annotation
    cluster_functions = annotate_gene_clusters(fasta_dir, index_file,
                                               cutoff=cutoff)

    # Read desired product annotation configuration file
    accept_list, ignore_list = fileio.read_functions_config_file(
                                                                functions_file)

    # Create a lookup map of clusters by function
    function_to_cluster_map = map_function_to_cluster(cluster_functions)

    # Isolate clusters whose consensus product annotation
    # meets the desired annotations from the config file
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
    # Prepare jobs for protein sequence alignment
    work_items = []
    for cluster in curated_clusters:
        cluster_path = fasta_dir.joinpath(cluster)
        outpath = output_dir.joinpath(cluster)

        records = [record for record in SeqIO.parse(cluster_path, "fasta")]

        if len(records) < min_size:
            continue

        work_items.append((records, outpath, cluster_functions[cluster]))

    # Run jobs
    multiprocess.parallelize(work_items, cores, dump_named_alignment,
                             verbose=verbose)


def dump_named_alignment(records, fasta_path, name):
    """Function to align a multiple sequence file and provide a name
    according to the specifications for hhsuite database building.

    :param records: SeqRecord protein sequence objects to be aligned
    :type records: list
    :param fasta_path: Path where the alignment file will be written
    :type fasta_path: pathlib.Path
    :param name: Name given to the seqeunce alignment
    :type name: str
    """
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


def annotate_gene_clusters(fasta_dir, index_file,
                           cutoff=DEFAULTS["cutoff"]):
    """Function to annotate a protein sequence cluster based on
    the product annotations of its members

    :param fasta_dir: Path to dir containing protein cluster fasta-alignments
    :type fasta_dir: pathlib.Path
    :param index_file: Path to gene index and relevant metadata file
    :type index_file: pathlib.Path
    :param cutoff: Frequency fraction required for the consensus product
    :type cutoff: float
    """
    gene_index = fileio.read_gene_index_file(index_file)

    cluster_functions = dict()
    for input_file in fasta_dir.iterdir():
        records = [record for record in SeqIO.parse(input_file, "fasta")]

        product_counts = {}
        for record in records:
            index = record.id

            gene_data = gene_index.get(index)
            if gene_data is None:
                continue

            product_count = product_counts.get(gene_data["product"], 0)
            product_count += 1
            product_counts[gene_data["product"]] = product_count

        products = [product for product, count in product_counts.items()
                    if (count / len(records)) >= cutoff]
        products.sort(key=lambda x: product_counts[x], reverse=True)

        cluster_product = HHSUITEDB_DEFAULTS["default_product"] 
        for product in products:
            if product == HHSUITEDB_DEFAULTS["default_product"]:
                continue

            cluster_product = product
            break

        cluster_functions[input_file.name] = cluster_product

    return cluster_functions


def main(unparsed_args):
    args = parse_curate_functions(unparsed_args)
    curate_gene_clusters(
                     args.input_dir, args.index_file, args.functions_config,
                     args.output_dir, cores=args.cpus, verbose=args.verbose,
                     accept_all=args.accept_all, min_size=args.min_size,
                     cutoff=args.cutoff)


if __name__ == "__main__":
    main(sys.argv[1:])
