"""Pipeline to compile the complete set of databases and data structures
required to run DEPhT from bacterial and phage sequences provided.
"""
import argparse
import json
import pathlib
import shutil
import sys

from Bio import SeqIO

from depht.__main__ import cleanup_flatfile_records, DEPHT_DIR
from depht.functions.annotation import annotate_record
from depht.functions.sniff_format import sniff_format
from depht_utils import PACKAGE_DIR
from depht_utils.data.defaults import (HHSUITEDB_DEFAULTS,
                                       MODEL_SCHEMA_DEFAULTS,
                                       SHELL_DB_DEFAULTS)
from depht_utils.functions import fileio
from depht_utils.pipelines.annotate_gene_clusters import annotate_gene_clusters
from depht_utils.pipelines.curate_functions import curate_functions
from depht_utils.pipelines.build_functions_db import build_functions_db
from depht_utils.pipelines.build_reference_db import build_reference_db
from depht_utils.pipelines.index_functions import index_functions
from depht_utils.pipelines.phamerate import execute_phamerate_pipeline
from depht_utils.pipelines.screen_conserved_phams import screen_conserved_phams
from depht_utils.pipelines.train_model import train_model


MODELS_DIR = DEPHT_DIR.joinpath("models")
DEFAULT_CONFIG = PACKAGE_DIR.joinpath("data/defaults.json")


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
    parser.add_argument("-c", "--config", type=pathlib.Path,
                        default=DEFAULT_CONFIG)
    parser.add_argument("-f", "--force", action="store_true")
    parser.add_argument("-a", "--auto_annotate", action="store_true")

    args = parser.parse_args()
    return args


def main(unparsed_args):
    """Main function for the command-line interface of the create_model
    pipeline.
    """
    args = parse_args(unparsed_args)

    create_model(args.model_name,
                 args.phage_sequences, args.bacterial_sequences,
                 verbose=args.verbose, config_file=args.config,
                 force=args.force, annotate=args.auto_annotate)


def load_config(config_file):
    """Function to load the master configuration file for the create_model
    pipeline and to verify its structure/contents.
    """
    with open(config_file, "r") as filehandle:
        config = json.load(filehandle)

    phameration_config = config.get("phameration")
    if phameration_config is None or not isinstance(phameration_config, dict):
        return None

    function_config = config.get("function")
    if function_config is None or not isinstance(function_config, dict):
        return None

    return config


def create_model(model_name, phage_sequences, bacterial_sequences,

                 verbose=False, config_file=DEFAULT_CONFIG, force=False,
                 annotate=False):
    # Load master configuration file, which contains paths for
    # configuration of sub-pipelines
    config = load_config(config_file)
    if config is None:
        print(f"Specified configuration file at {config_file} "
              "is incorrectly formatted.\n"
              "Please check the formatting of this file before continuing.")
        return

    # Create local model directory
    model_dir = MODELS_DIR.joinpath(model_name)
    if model_dir.is_dir():
        if force:
            shutil.rmtree(model_dir)
        else:
            print("There already exists a DEPHT model with the name "
                  f"'{model_name}'.\n"
                  "Please rename your model or use the -f flag to forcibly "
                  "overwrite the already existing model.")
            return
    model_dir.mkdir()

    # Create temporary working directory
    tmp_dir = model_dir.joinpath("tmp")
    tmp_dir.mkdir()

    if verbose:
        print("Collecting/annotating bacterial sequences...")
    bacterial_tmp = tmp_dir.joinpath("bacteria")
    bacterial_tmp.mkdir()
    bacterial_fasta_files, bacterial_gb_files, num_bacteria = clean_sequences(
                                          bacterial_sequences, bacterial_tmp,
                                          annotate=annotate)

    if num_bacteria == 0:
        print("Could not recognize any bacterial sequence files "
              "in the specified directory.\n"
              "Please check the formatting and contents of your files.")
        shutil.rmtree(model_dir)
        return

    if verbose:
        print("Collecting/annotating phage sequences...")
    phage_tmp = tmp_dir.joinpath("phage")
    phage_tmp.mkdir()
    phage_fasta_files, phage_gb_files, num_phages = clean_sequences(
                                                phage_sequences, phage_tmp,
                                                annotate=annotate)

    if num_phages == 0:
        print("Could not recognize any phage sequence files "
              "in the specified directory.\n"
              "Please check the formatting and contents of your files.")
        shutil.rmtree(model_dir)
        return

    if verbose:
        print("Training phage/bacterial classifier...")
    train_model(model_name, phage_fasta_files, bacterial_fasta_files)

    if verbose:
        print("Building bacterial reference database...")
    reference_db_dir = model_dir.joinpath(
                                    MODEL_SCHEMA_DEFAULTS["reference_db"])
    reference_db_dir.mkdir()
    build_reference_db(bacterial_fasta_files, reference_db_dir)

    if verbose:
        print("Building shell genome content database...")
    shell_db_tmp = tmp_dir.joinpath(MODEL_SCHEMA_DEFAULTS["shell_db"])
    shell_db_tmp.mkdir()
    shell_db_dir = model_dir.joinpath(
                                    MODEL_SCHEMA_DEFAULTS["shell_db"])
    shell_db_dir.mkdir()
    create_shell_db(bacterial_gb_files, shell_db_dir, config, shell_db_tmp,
                    verbose=verbose)

    if verbose:
        print("Building phage homolog profile database...")
    phage_homologs_tmp = tmp_dir.joinpath(
                                    MODEL_SCHEMA_DEFAULTS["phage_homologs_db"])
    phage_homologs_tmp.mkdir()
    phage_homologs_dir = model_dir.joinpath(
                                    MODEL_SCHEMA_DEFAULTS["phage_homologs_db"])
    phage_homologs_dir.mkdir()
    create_phage_homologs_db(phage_sequences, phage_homologs_dir, config,
                             phage_homologs_tmp, verbose=verbose)

    shutil.rmtree(tmp_dir)


def clean_sequences(input_dir, output_dir, annotate=False):
    fasta_dir = output_dir.joinpath("fasta")
    fasta_dir.mkdir() 

    gb_dir = output_dir.joinpath("gb")
    gb_dir.mkdir()

    seq_count = 0
    for input_file in input_dir.iterdir():
        file_fmt = sniff_format(input_file) 
        if file_fmt is None:
            continue

        records = [record for record in SeqIO.parse(input_file, file_fmt)]
   
        if file_fmt == "fasta" or annotate:
            for record in records:
                annotate_record(record, output_dir)
        else:
            cleanup_flatfile_records(records)

        fasta_file = fasta_dir.joinpath(".".join([input_file.stem, "fasta"]))
        SeqIO.write(records, fasta_file, "fasta")

        gb_file = gb_dir.joinpath(".".join([input_file.stem, "gb"]))
        SeqIO.write(records, gb_file, "gb")

        seq_count += 1

    return (fasta_dir, gb_dir, seq_count)


def create_shell_db(bacterial_sequences, output_dir, config, tmp_dir,
                    verbose=False):
    if verbose:
        print("...indexing bacterial protein sequences...")
    # Create a simple fasta-based database from the given bacterial sequences
    fasta_file, index_file, cluster_file = index_functions(
                                                bacterial_sequences, tmp_dir,
                                                name=SHELL_DB_DEFAULTS["name"])

    gene_clusters_dir = tmp_dir.joinpath("phams")
    gene_clusters_dir.mkdir()
    # Phamerate bacterial sequences

    if verbose:
        print("...clustering bacterial protein sequences...")
    phamerate_config = PACKAGE_DIR.joinpath(config["phameration"]["bacteria"])
    execute_phamerate_pipeline(fasta_file, gene_clusters_dir, phamerate_config)

    if cluster_file is None:
        cluster_file = create_cluster_schema(index_file, gene_clusters_dir,
                                             tmp_dir)

    if verbose:
        print("...screening for shell genome content...")
    screen_conserved_phams(gene_clusters_dir, output_dir, index_file,
                           cluster_file)

    fasta_file.replace(output_dir.joinpath(fasta_file.name))


def create_phage_homologs_db(phage_sequences, output_dir, config, tmp_dir,
                             verbose=False):
    if verbose:
        print("...indexing phage protein sequences...")
    # Create a simple fasta-based database from the given phage sequences
    fasta_file, index_file, cluster_file = index_functions(
                                                phage_sequences, tmp_dir)

    gene_clusters_dir = tmp_dir.joinpath("phams")
    gene_clusters_dir.mkdir()

    if verbose:
        print("...clustering phage protein sequences...")
    # Phamerate bacterial sequences
    phamerate_config = PACKAGE_DIR.joinpath(config["phameration"]["bacteria"])
    execute_phamerate_pipeline(fasta_file, gene_clusters_dir, phamerate_config)

    if verbose:
        print("...assigning annotations to phage protein clusters...")
    annotation_index = annotate_gene_clusters(gene_clusters_dir, index_file,
                                              tmp_dir)

    if verbose:
        print("...curating for essential phage protein clusters...")
    curated_clusters_dir = tmp_dir.joinpath("curated_phams")
    curated_clusters_dir.mkdir()
    functions_config = PACKAGE_DIR.joinpath(config["functions"]["essential"])
    curate_functions(gene_clusters_dir, annotation_index, functions_config,
                     curated_clusters_dir, verbose=verbose, cores=6)

    if verbose:
        print("...creating database of essential phage protein profiles...")
    build_functions_db(curated_clusters_dir, output_dir, name="essential")

    if verbose:
        print("...curating for accessory phage protein clusters...")
    extended_curated_clusters_dir = tmp_dir.joinpath("extended_curated_phams")
    extended_curated_clusters_dir.mkdir()
    functions_config = PACKAGE_DIR.joinpath(config["functions"]["extended"])
    curate_functions(gene_clusters_dir, annotation_index, functions_config,
                     curated_clusters_dir, verbose=verbose, accept_all=True,
                     cores=6)

    if verbose:
        print("...creating database of accessory phage protein profiles...")
    build_functions_db(curated_clusters_dir, output_dir, name="extended")


def create_cluster_schema(index_file, gene_clusters_dir, tmp_dir):
    gene_index = fileio.read_gene_index_file(index_file)

    clustered_ids = list()
    ids = list()
    for data_dict in list(gene_index.values()):
        ids.append(data_dict["parent"])

    clustered_ids.append(ids)

    cluster_file = tmp_dir.joinpath(".".join([HHSUITEDB_DEFAULTS["name"],
                                              "ci"])) 
    fileio.write_cluster_file(clustered_ids, cluster_file)

    return cluster_file

    
    pass



if __name__ == "__main__":
    unparsed_args = sys.argv
    if len(unparsed_args) <= 1:
        unparsed_args.append("-h")

    main(unparsed_args[1:])
