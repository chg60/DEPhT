"""Pipeline to compile the complete set of databases and data structures
required to run DEPhT from bacterial and phage sequences provided.
"""
import argparse
import json
import pathlib
import shutil
import sys

from Bio import SeqIO
from phamclust.__main__ import get_clusters

from depht.__main__ import MODEL_DIR
from depht.data import GLOBAL_VARIABLES
from depht.functions.annotation import (annotate_record,
                                        cleanup_flatfile_records)
from depht.functions.fasta import parse_fasta
from depht.functions.multiprocess import CPUS, parallelize
from depht.functions.sniff_format import sniff_format
from depht_train.data import PARAMETERS
from depht_train.functions import fileio
from depht_train.pipelines.build_HMM_db import build_HMM_db
from depht_train.pipelines.build_reference_db import build_reference_db
from depht_train.pipelines.curate_gene_clusters import (
    curate_gene_clusters)
from depht_train.pipelines.index_sequences import index_sequences
from depht_train.pipelines.phamerate import (execute_phamerate_pipeline,
                                             parse_param_dict)
from depht_train.pipelines.screen_conserved_phams import (
    screen_conserved_phams)
from depht_train.pipelines.train_model import train_model


def parse_args(unparsed_args):
    """Function to parse command line arguments for running the create_model
    pipeline.

    :param unparsed_args: List of command line arguments
    :type unparsed_args: list
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("name", type=str,
                        help="a name for the newly created model")
    parser.add_argument("phages", type=pathlib.Path,
                        help="directory containing phage sequences to train "
                             "with")
    parser.add_argument("bacteria", type=pathlib.Path,
                        help="directory containing bacterial sequences to "
                             "train with")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="print progress messages to console")
    parser.add_argument("-c", "--cpus", type=int, default=CPUS, metavar="",
                        help=f"number of CPU cores to use [default: {CPUS}]")
    parser.add_argument("-f", "--force", action="store_true",
                        help="overwrite existing model with the same name")

    parser.add_argument("--config", type=pathlib.Path, metavar="",
                        help=f"configuration file to customize pipeline "
                             f"parameters")
    parser.add_argument("--bacteria-clusters", type=pathlib.Path, metavar="",
                        help="a CSV file with manually determined clusters "
                             "for the input bacterial sequences")
    parser.add_argument("--prophage-coords", type=pathlib.Path, metavar="",
                        help="a CSV file with the coordinates of any known "
                             "prophages in the input bacterial sequences")

    return parser.parse_args(unparsed_args)


def load_config(config_file):
    """Function to load the master configuration file for the create_model
    pipeline and to verify its structure/contents.
    """
    config = PARAMETERS
    if config_file is None:
        return config
    elif not config_file.is_file():
        return config

    with open(config_file, "r") as filehandle:
        input_config = json.load(filehandle)

    config.update(input_config)
    return config


def create_model(name, phages, bacteria, bact_clusters=None,
                 prophage_coords=None, config=None, force=False,
                 verbose=False, cpus=1):
    """Function to create a DEPhT model from bacterial and phage sequences.

    :param name: Name given to the created model.
    :type name: str
    :param phages: Directory of phage sequence files.
    :type phages: pathlib.Path
    :param bacteria: Directory of bacterial sequence files.
    :type bacteria: pathlib.Path
    :param verbose: Boolean to toggle pipeline verbosity.
    :type verbose: bool
    :param config: Config file to overwrite pipeline parameters.
    :type config: pathlib.Path
    :param force: Allows pipeline to overwrite already existing models:
    :type force: bool
    :param cpus: Number of CPUs to use during the pipeline.
    :type cpus: int
    :param bact_clusters: File describing bacteria clade membership
    :type bact_clusters: pathlib.Path
    :param prophage_coords: File describing known prophage content in bacteria
    :type prophage_coords: pathlib.Path
    """
    # Load master configuration file, which contains paths for
    # configuration of sub-pipelines
    config = load_config(config)
    if config is None:
        print(f"Specified configuration file at {config} "
              "is incorrectly formatted or contains invalid information.\n"
              "Please check the formatting/contents of this file before "
              "continuing.")
        return

    # Proactively create local model directory structure
    dir_map = create_model_structure(name, force=force)
    if dir_map is None:
        print("There already exists a DEPHT model with the name "
              f"'{name}'.\n"
              "Please rename your model or use the -f flag to forcibly "
              "overwrite the already existing model.")
        return

    # Collect bacterial sequences and write fasta and genbank files
    if verbose:
        print("Collecting/annotating bacterial sequences...")
    bacterial_data_tuple = collect_sequences(
        bacteria, dir_map["bacterial_tmp"],
        GLOBAL_VARIABLES["bacterial_sequences"]["name"],
        config["bacterial_sequences"],
        annotate=True, verbose=verbose, cpus=cpus,
        cluster_table=bact_clusters)

    if bacterial_data_tuple[2] == 0:
        print("Could not recognize any bacterial sequence files "
              "in the specified directory.\n"
              "Please check the formatting and contents of your files.")
        shutil.rmtree(dir_map["model_dir"])
        return

    # Collect phage sequences and write fasta and genbank files
    if verbose:
        print("Collecting/annotating phage sequences...")
    phage_data_tuple = collect_sequences(
        phages, dir_map["phage_tmp"],
        GLOBAL_VARIABLES["phage_sequences"]["name"],
        config["phage_sequences"],
        verbose=verbose, cpus=cpus, cluster=False,
        singletons=True)

    if phage_data_tuple[2] == 0:
        print("Could not recognize any phage sequence files "
              "in the specified directory.\n"
              "Please check the formatting and contents of your files.")
        shutil.rmtree(dir_map["model_dir"])
        return

    if verbose:
        print("Training phage/bacterial classifier...")
    train_model(name, phage_data_tuple[1], bacterial_data_tuple[1],
                cpus=cpus, prophages=prophage_coords)

    if verbose:
        print("\nBuilding bacterial reference database...")
    build_reference_db(bacterial_data_tuple[0], dir_map["reference_db_dir"])

    if verbose:
        print("\nBuilding shell genome content database...")
    rep_threshold = config["shell_db"]["rep_threshold"]
    screen_conserved_phams(bacterial_data_tuple[5], dir_map["shell_db_dir"],
                           bacterial_data_tuple[4], bacterial_data_tuple[6],
                           rep_threshold=rep_threshold)

    # Moves bacterial gene product fasta file into the shell genome database
    bacterial_data_tuple[3].replace(dir_map["shell_db_dir"].joinpath(
                                                bacterial_data_tuple[3].name))

    if verbose:
        print("\nBuilding phage homolog profile database...")
    create_phage_homologs_db(phage_data_tuple, dir_map["phage_homologs_dir"],
                             config, dir_map["phage_homologs_tmp"],
                             verbose=verbose, cpus=cpus)

    # Cleans up after the pipeline
    shutil.rmtree(dir_map["tmp_dir"])


def create_model_structure(model_name, force=False):
    """Function to create a DEPhT model's directory structure.

    :param model_name: Name given to the created model.
    :type model_name: str
    :param force: Allows pipeline to overwrite already existing models.
    :type force: bool
    """
    dir_map = dict()

    # Create local model directory
    model_dir = MODEL_DIR.joinpath(model_name)
    dir_map["model_dir"] = model_dir
    if model_dir.is_dir():
        if force:
            shutil.rmtree(model_dir)
        else:
            return None

    model_dir.mkdir()

    # Create temporary working directory
    tmp_dir = model_dir.joinpath("tmp")
    dir_map["tmp_dir"] = tmp_dir
    tmp_dir.mkdir()

    # Create temporary bacterial sequence directory
    bacterial_tmp = tmp_dir.joinpath("bacteria")
    dir_map["bacterial_tmp"] = bacterial_tmp
    bacterial_tmp.mkdir()

    # Create temporary phage sequence directory
    phage_tmp = tmp_dir.joinpath("phage")
    dir_map["phage_tmp"] = phage_tmp
    phage_tmp.mkdir()

    # Create model reference blast database directory
    reference_db_dir = model_dir.joinpath(
                                GLOBAL_VARIABLES["reference_db"]["dir_name"])
    dir_map["reference_db_dir"] = reference_db_dir
    reference_db_dir.mkdir()

    # Create model shell database directory and temporary directory
    shell_db_tmp = tmp_dir.joinpath(
                                GLOBAL_VARIABLES["shell_db"]["dir_name"])
    dir_map["shell_db_tmp"] = shell_db_tmp
    shell_db_tmp.mkdir()
    shell_db_dir = model_dir.joinpath(
                                GLOBAL_VARIABLES["shell_db"]["dir_name"])
    dir_map["shell_db_dir"] = shell_db_dir
    shell_db_dir.mkdir()

    # Create model phage homolog HMM profile directory and temporary directory
    phage_homologs_tmp = tmp_dir.joinpath(
                                GLOBAL_VARIABLES["phage_homologs"]["dir_name"])
    dir_map["phage_homologs_tmp"] = phage_homologs_tmp
    phage_homologs_tmp.mkdir()
    phage_homologs_dir = model_dir.joinpath(
                                GLOBAL_VARIABLES["phage_homologs"]["dir_name"])
    dir_map["phage_homologs_dir"] = phage_homologs_dir
    phage_homologs_dir.mkdir()
    return dir_map


def collect_sequences(sequences, output_dir, name, config,
                      annotate=False, verbose=False, cpus=1,
                      singletons=False, cluster=True, cluster_table=None):
    """Function to annotate, phamerate gene products, and cluster input
    sequences.

    :param sequences: Directory of sequence files.
    :type sequences: pathlib.Path
    :param output_dir: Directory to write sequence files to.
    :type output_dir: pathlib.Path
    :param name: Name to give to the stem of the created sequence files.
    :type name: str
    :param config: Config dictionary that contains pipeline parameters
    :type config: dict
    :param annotate: Controls whether all sequence files are auto-annotated
    :type annotate: bool
    :param verbose: Boolean to toggle pipeline verbosity.
    :type verbose: bool
    :param cpus:
    :type cpus: int
    :param singletons: Controls whether singleton clusters are retained.
    :type singletons: bool
    :param cluster: Controls whether sequences are clustered.
    :type cluster: bool
    :param cluster_table: File describing sequence clade membership
    :type cluster_table: pathlib.Path
    """
    fasta_sequences, gb_sequences, seq_count = clean_sequences(
                                    sequences, output_dir, annotate=annotate,
                                    verbose=verbose, cpus=cpus)

    fasta_file, index_file, cluster_file = index_sequences(
                                    sequences, output_dir,
                                    name=name, cluster_table=cluster_table)

    gene_clusters_dir = output_dir.joinpath("phams")

    if verbose:
        print("...clustering gene product sequences...")
    phamerate_config = config["phameration"]
    first_iter_params, second_iter_params = parse_param_dict(phamerate_config)
    execute_phamerate_pipeline(fasta_file, gene_clusters_dir,
                               first_iter_params, second_iter_params)

    clusters = None
    if cluster_file is None and cluster:
        if verbose:
            print("...clustering genomes based on shared gene content...")
        cluster_file, clusters = create_cluster_schema(
                                            index_file, gene_clusters_dir,
                                            output_dir, name, config,
                                            singletons=singletons)

    return (fasta_sequences, gb_sequences, seq_count, fasta_file,
            index_file, gene_clusters_dir, cluster_file, clusters)


def clean_sequences(input_dir, output_dir, annotate=False, verbose=False,
                    cpus=CPUS, trna=False):
    """Formats sequence files, and annotates them as applicable.

    :param input_dir: Directory of sequence files
    :type input_dir: pathlib.Path
    :param output_dir: Directory to write sequence files to.
    :type output_dir: pathlib.Path
    :param annotate: Controls whether sequence files are auto-annotated.
    :type annotate: bool
    :param verbose: Boolean to toggle pipeline verbosity.
    :type verbose: bool
    :param cpus: Number of CPUs to use during the pipeline.
    :type cpus: int
    :param trna: Controls whether sequence files are annotated with tRNAs.
    :type trna: bool
    """
    fasta_dir = output_dir.joinpath("fasta")
    fasta_dir.mkdir()

    gb_dir = output_dir.joinpath("gb")
    gb_dir.mkdir()

    work_items = list()
    for input_file in input_dir.iterdir():
        work_items.append((input_file, output_dir, fasta_dir, gb_dir,
                           annotate, trna))

    seq_count = sum(parallelize(work_items, cpus, clean_sequence,
                                verbose=verbose))

    return fasta_dir, gb_dir, seq_count


def clean_sequence(input_file, output_dir, fasta_dir, gb_dir, annotate=False,
                   trna=False):
    """Process function to format a sequence file and annotates it as
    applicable.

    :param input_file: Directory of sequence files
    :type input_file: pathlib.Path
    :param output_dir:
    :type output_dir:
    :param fasta_dir: Directory to write the fasta-formatted sequence file to.
    :type fasta_dir: pathlib.Path
    :param gb_dir: Directory to write the gb-formatted sequence file to.
    :type gb_dir: pathlib.Path
    :param annotate: Controls whether sequence files are auto-annotated.
    :type annotate: bool
    :param trna: Controls whether sequence files are annotated with tRNAs.
    :type trna: bool
    """
    file_fmt = sniff_format(input_file)
    if file_fmt not in ("fasta", "genbank"):
        return 0

    records = [record for record in SeqIO.parse(input_file, file_fmt)]

    if file_fmt == "fasta" or annotate:
        for record in records:
            record.features = list()
            annotate_record(record, output_dir, trna=trna)

    cleanup_flatfile_records(records)

    fasta_file = fasta_dir.joinpath(".".join([input_file.stem, "fasta"]))
    SeqIO.write(records, fasta_file, "fasta")

    gb_file = gb_dir.joinpath(".".join([input_file.stem, "gb"]))
    SeqIO.write(records, gb_file, "gb")

    return 1


def create_phage_homologs_db(phage_data_tuple, output_dir, config, tmp_dir,
                             cpus=CPUS, verbose=False):
    """

    :param phage_data_tuple:
    :param output_dir:
    :param config:
    :param tmp_dir:
    :param cpus:
    :param verbose:
    :return:
    """
    min_hmm_count = config["phage_homologs"]["min_hmm_count"]

    if verbose:
        print("...curating for essential phage protein clusters...")
    curated_clusters_dir = tmp_dir.joinpath("curated_phams")
    curated_clusters_dir.mkdir()
    accepted_functions = config["phage_homologs"][
                                "essential_annotations"]["LIKE"]
    ignored_functions = config["phage_homologs"][
                                "essential_annotations"]["NOT LIKE"]
    curate_gene_clusters(phage_data_tuple[5], phage_data_tuple[4],
                         curated_clusters_dir,
                         accepted_functions, ignored_functions,
                         verbose=verbose, cores=cpus,
                         min_hmm_count=min_hmm_count)

    if verbose:
        print("...creating database of essential phage protein profiles...")
    build_HMM_db(curated_clusters_dir, output_dir, name="essential")

    if verbose:
        print("...curating for accessory phage protein clusters...")
    extended_curated_clusters_dir = tmp_dir.joinpath("extended_curated_phams")
    extended_curated_clusters_dir.mkdir()
    accepted_functions = config["phage_homologs"][
                                "extended_annotations"]["LIKE"]
    ignored_functions = config["phage_homologs"][
                               "extended_annotations"]["NOT LIKE"]
    curate_gene_clusters(phage_data_tuple[5], phage_data_tuple[4],
                         curated_clusters_dir,
                         accepted_functions, ignored_functions,
                         verbose=verbose,
                         accept_all=True, cores=cpus,
                         min_hmm_count=min_hmm_count)

    if verbose:
        print("...creating database of accessory phage protein profiles...")
    build_HMM_db(curated_clusters_dir, output_dir, name="extended")


def create_cluster_schema(index_file, gene_clusters_dir, tmp_dir, name,
                          config, singletons=False):
    """Create the cluster schema for a group of sequences.

    :param index_file: Tsv table index of sequence gene products
    :type index_file: pathlib.Path
    :param gene_clusters_dir: Directory of pham (clustered gene product) fastas
    :type gene_clusters_dir: pathlib.Path
    :param tmp_dir: Directory to put pipeline temporary data.
    :type  tmp_dir: pathlib.Path
    :param name: Name to give to the stem of the created sequence files.
    :type name: str
    :param config: Config  dictionary that contains pipeline parameters
    :type config: dict
    :param singletons: Controls whether singleton clusters are retained.
    :type singletons: bool
    """
    gene_index = fileio.read_gene_index_file(index_file)

    genome_pham_tuples = list()
    for pham_file in gene_clusters_dir.iterdir():
        if pham_file.suffix != ".fasta":
            continue

        pham = pham_file.stem
        pham_headers, pham_sequences = parse_fasta(pham_file)

        for header in pham_headers:
            parent_id = gene_index[header]["parent"]
            genome_pham_tuples.append([parent_id, pham])

    genome_pham_tsv_file = tmp_dir.joinpath(".".join([name, "tsv"]))
    with open(genome_pham_tsv_file, "w") as filehandle:
        for genome_pham_tuple in genome_pham_tuples:
            line = "\t".join(genome_pham_tuple)
            filehandle.write(f"{line}\n")

    clusters = get_clusters(genome_pham_tsv_file,
                            config["pham_clust"]["gcs_threshold"])

    clustered_ids = list()
    for cluster_matrix in clusters:
        if len(cluster_matrix.node_names) <= 1 and not singletons:
            continue

        clustered_ids.append(cluster_matrix.node_names)

    cluster_file = tmp_dir.joinpath(".".join([name, "ci"]))
    fileio.write_cluster_file(clustered_ids, cluster_file)

    return cluster_file, clusters


def main(unparsed_args=None):
    """Main function for the command-line interface of the create_model
    pipeline.
    """
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])

    create_model(args.name, args.phages, args.bacteria,
                 bact_clusters=args.bacteria_clusters,
                 prophage_coords=args.prophage_coords,
                 config=args.config, verbose=args.verbose,
                 force=args.force, cpus=args.cpus)


if __name__ == "__main__":
    main()
