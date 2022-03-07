"""Pipeline to compile the complete set of databases and data structures
required to run DEPhT from bacterial and phage sequences provided.
"""
import argparse
import json
import pathlib
import shutil

from Bio import SeqIO
from phamclust.__main__ import get_clusters

from depht.__main__ import MODEL_DIR
from depht.data import GLOBAL_VARIABLES
from depht.functions.annotation import (annotate_record,
                                        cleanup_flatfile_records)
from depht.functions.fasta import parse_fasta
from depht.functions.multiprocess import CPUS, parallelize
from depht.functions.sniff_format import sniff_format
from depht_utils.data import PARAMETERS
from depht_utils.functions import fileio
from depht_utils.pipelines.build_HMM_db import build_HMM_db
from depht_utils.pipelines.build_reference_db import build_reference_db
from depht_utils.pipelines.curate_gene_clusters import (
    curate_gene_clusters)
from depht_utils.pipelines.index_sequences import index_sequences
from depht_utils.pipelines.phamerate import (execute_phamerate_pipeline,
                                             parse_param_dict)
from depht_utils.pipelines.screen_conserved_phams import (
    screen_conserved_phams)
from depht_utils.pipelines.train_model import train_model


def parse_args(unparsed_args):
    """Function to parse command line arguments for running the create_model
    pipeline.

    :param unparsed_args: List of command line arguments
    :type unparsed_args: list
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", "--model_name", type=str, required=True)
    parser.add_argument("-p", "--phage_sequences", type=pathlib.Path,
                        required=True)
    parser.add_argument("-b", "--bacterial_sequences", type=pathlib.Path,
                        required=True)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--cpus", type=int, default=CPUS)
    parser.add_argument("-c", "--config", type=pathlib.Path,
                        default=None)
    parser.add_argument("-f", "--force", action="store_true")
    parser.add_argument("-a", "--auto_annotate", action="store_true")

    args = parser.parse_args(unparsed_args)
    return args


def main(unparsed_args):
    """Main function for the command-line interface of the create_model
    pipeline.
    """
    args = parse_args(unparsed_args)

    create_model(args.model_name,
                 args.phage_sequences, args.bacterial_sequences,
                 verbose=args.verbose, config_file=args.config,
                 force=args.force, annotate=args.auto_annotate, cpus=args.cpus)


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


def create_model(model_name, phage_sequences, bacterial_sequences,
                 verbose=False, config_file=None, force=False,
                 cpus=1, annotate=False):
    # Load master configuration file, which contains paths for
    # configuration of sub-pipelines
    config = load_config(config_file)
    if config is None:
        print(f"Specified configuration file at {config_file} "
              "is incorrectly formatted or contains invalid information.\n"
              "Please check the formatting/contents of this file before "
              "continuing.")
        return

    # Proactively create local model directory structure
    dir_map = create_model_structure(model_name, force=force)
    if dir_map is None:
        print("There already exists a DEPHT model with the name "
              f"'{model_name}'.\n"
              "Please rename your model or use the -f flag to forcibly "
              "overwrite the already existing model.")
        return

    # Annotate bacterial sequences and write fasta and genbank files
    if verbose:
        print("Collecting/annotating bacterial sequences...")
    bacterial_data_tuple = collect_sequences(
                            bacterial_sequences, dir_map["bacterial_tmp"],
                            GLOBAL_VARIABLES["bacterial_sequences"]["name"],
                            config["bacterial_sequences"],
                            annotate=False, verbose=verbose, cpus=cpus)

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
                            phage_sequences, dir_map["phage_tmp"],
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
    train_model(model_name, phage_data_tuple[1], bacterial_data_tuple[1],
                cpus=cpus)

    if verbose:
        print("\nBuilding bacterial reference database...")
    build_reference_db(bacterial_data_tuple[0], dir_map["reference_db_dir"])

    if verbose:
        print("\nBuilding shell genome content database...")
    rep_threshold = config["shell_db"]["rep_threshold"]
    screen_conserved_phams(bacterial_data_tuple[5], dir_map["shell_db_dir"],
                           bacterial_data_tuple[4], bacterial_data_tuple[6],
                           rep_threshold=rep_threshold)

    bacterial_data_tuple[3].replace(dir_map["shell_db_dir"].joinpath(
                                                bacterial_data_tuple[3].name))

    if verbose:
        print("\nBuilding phage homolog profile database...")
    create_phage_homologs_db(phage_data_tuple, dir_map["phage_homologs_dir"],
                             config, dir_map["phage_homologs_tmp"],
                             verbose=verbose, cpus=cpus)

    shutil.rmtree(dir_map["tmp_dir"])


def create_model_structure(model_name, force=False):
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
                      annotate=False, verbose=False, cpus=CPUS,
                      singletons=False, cluster=True):
    fasta_sequences, gb_sequences, seq_count = clean_sequences(
                                    sequences, output_dir, annotate=annotate,
                                    verbose=verbose, cpus=cpus)

    fasta_file, index_file, cluster_file = index_sequences(
                                    sequences, output_dir,
                                    name=name)

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

    return (fasta_dir, gb_dir, seq_count)


def clean_sequence(input_file, output_dir, fasta_dir, gb_dir, annotate=False,
                   trna=False):
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
