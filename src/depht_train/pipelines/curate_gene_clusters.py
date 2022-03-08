"""Curate phage gene phamilies, looking for conserved functional
annotations.
"""
import argparse
import pathlib
import sys

from Bio import SeqIO, AlignIO

from depht.data import GLOBAL_VARIABLES
from depht.functions.multiprocess import CPUS, parallelize
from depht_train.data import PARAMETERS
from depht_train.functions import clustalo, fileio

MIN_HMM_COUNT = PARAMETERS["phage_homologs"]["min_hmm_count"]
AC_THRESHOLD = PARAMETERS["phage_homologs"]["annotation_consensus_threshold"]


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_args(unparsed_args):
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)
    parser.add_argument("index_file", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--name", type=str,
                        default=GLOBAL_VARIABLES["phage_sequences"]["name"])
    parser.add_argument("-np", "--cpus", type=int, default=CPUS)
    parser.add_argument("-c", "--config", type=pathlib.Path, default=None)

    parser.add_argument("-ct", "--consensus_threshold", type=float,
                        default=AC_THRESHOLD)
    parser.add_argument("-all", "--accept_all", action="store_true")
    parser.add_argument("-mc", "--min_hmm_count", type=int)

    args = parser.parse_args(unparsed_args)
    return args


def curate_gene_clusters(fasta_dir, index_file, output_dir,
                         accepted_functions, ignored_functions,
                         accept_all=False, verbose=False, cores=1,
                         min_hmm_count=MIN_HMM_COUNT,
                         consensus_threshold=AC_THRESHOLD):
    """Annotate and curate gene clusters, discarding those gene clusters
    with undesired or unclear consensus functions

    :param accept_all:
    :param fasta_dir: Path to dir containing protein cluster fasta-alignments
    :type fasta_dir: pathlib.Path
    :param index_file: Path to gene index and relevant metadata file
    :type index_file: pathlib.Path
    :param accepted_functions: List of function annotations to accept
    :type accepted_functions: list[str]
    :param ignored_functions: List of function annotations to ignore
    :type ignored_functions: list[str]
    :param output_dir: Path to directory to dump named/aligned sequence files
    :type output_dir: pathlib.Path
    :param cores: Number of processes to utilize
    :type cores: int
    :param verbose: Toggle pipeline print statements
    :type verbose: bool
    :param min_hmm_count: Minimum size threshold for protein clusters
    :type min_hmm_count: int
    :param consensus_threshold: Fraction required for the consensus product.
    :type consensus_threshold: float
    """
    # Annotate gene clusters based on consensus product annotation
    cluster_functions = annotate_gene_clusters(
                                    fasta_dir, index_file,
                                    consensus_threshold=consensus_threshold)

    # Create a lookup map of clusters by function
    function_to_cluster_map = map_function_to_cluster(cluster_functions)

    # Isolate clusters whose consensus product annotation
    # meets the desired annotations from the config file
    curated_clusters = []
    for cluster_function, clusters in function_to_cluster_map.items():
        ignore = False
        for ignored_function in ignored_functions:
            if ignored_function.lower() in cluster_function.lower():
                ignore = True
                break

        if ignore:
            continue

        if not accept_all:
            for accepted_function in accepted_functions:
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

        if len(records) < min_hmm_count:
            continue

        work_items.append((records, outpath, cluster_functions[cluster]))

    # Run jobs
    parallelize(work_items, cores, dump_named_alignment, verbose=verbose)


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
    """

    :param cluster_functions:
    :return:
    """
    function_to_cluster_map = {}
    for cluster_name, cluster_function in cluster_functions.items():
        mapped_clusters = function_to_cluster_map.get(cluster_function, list())
        mapped_clusters.append(cluster_name)
        function_to_cluster_map[cluster_function] = mapped_clusters

    return function_to_cluster_map


def annotate_gene_clusters(fasta_dir, index_file,
                           consensus_threshold=AC_THRESHOLD):
    """Function to annotate a protein sequence cluster based on
    the product annotations of its members

    :param fasta_dir: Path to dir containing protein cluster fasta-alignments
    :type fasta_dir: pathlib.Path
    :param index_file: Path to gene index and relevant metadata file
    :type index_file: pathlib.Path
    :param consensus_threshold:  Fraction required for the consensus product
    :type consensus_threshold: float
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
                    if (count / len(records)) >= consensus_threshold]
        products.sort(key=lambda x: product_counts[x], reverse=True)

        cluster_product = GLOBAL_VARIABLES["sequences"]["default_product"]
        for product in products:
            if product == GLOBAL_VARIABLES["sequences"]["default_product"]:
                continue

            cluster_product = product
            break

        cluster_functions[input_file.name] = cluster_product

    return cluster_functions


def main(unparsed_args=None):
    """Commandline entrypoint to this module."""
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])

    accept_list = None
    ignore_list = None
    if args.config is not None:
        if args.config.is_file():
            # Read desired product annotation configuration file
            accept_list, ignore_list = fileio.read_functions_config_file(
                                                        args.config)

    if accept_list is None or ignore_list is None:
        if accept_list is None:
            accept_list = PARAMETERS["phage_homologs"][
                                     "essential_annotations"]["LIKE"]
        if ignore_list is None:
            ignore_list = PARAMETERS["phage_homologs"][
                                     "essential_annotations"]["NOT LIKE"]

    curate_gene_clusters(
        args.input_dir, args.index_file, args.output_dir,
        accept_list, ignore_list, cores=args.cpus, verbose=args.verbose,
        accept_all=args.accept_all, min_hmm_count=args.min_hmm_count,
        consensus_threshold=args.consensus_threshold)


if __name__ == "__main__":
    main()
