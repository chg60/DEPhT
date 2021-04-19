"""
Module for training a Random Forest Classifier with best accuracy
in discriminating between bacterial and prophage genes.

For best results, limit the scope to a single bacterial genus or
phylum.
"""

import sys
import argparse
import pathlib
import pickle

import pandas as pd
from Bio import SeqIO

from prophicient.__main__ import MIN_LENGTH, META_LENGTH
from prophicient.functions.gene_prediction import predict_cds_features
from prophicient.functions.prophage_prediction import contig_to_dataframe, MODEL_PATH
from train_model.functions.model_training import train_random_forest_classifier

EPILOG = """
Bacterial genomes can be incomplete assemblies, so long as there are no
ambiguous nucleotides. In this case, training will be done on all contigs
longer than 20,000 bp. Bacterial genomes used for training should not have
any prophages.
"""


def parse_args(arguments):
    """
    Uses argparse module to parse the command line arguments.

    :param arguments: command line arguments
    :type arguments: list
    :return: parsed_args
    """
    p = argparse.ArgumentParser(description=__doc__, epilog=EPILOG)
    p.add_argument("prophages", type=pathlib.Path,
                   help="path to a directory containing prophage FASTA files")
    p.add_argument("bacteria", type=pathlib.Path,
                   help="path to a directory containing bacterial FASTA files")
    p.add_argument("--verbose", action="store_true",
                   help="print progress messages as program runs")
    return p.parse_args(arguments)


def annotate_sequences(directory):
    """
    Convenience method to facilitate annotation of many sequences
    in a directory.

    :param directory: path where FASTA files live
    :type directory: pathlib.Path
    :return: dataframe
    """
    dfs = list()

    for fasta_path in directory.iterdir():
        # Check file extension to see if this file is a FASTA
        if fasta_path.suffix not in (".fasta", ".fna"):
            # print(f"'{str(fasta_path)}' appears not to be a FASTA file...")
            continue

        # Parse any files that have a valid FASTA extension
        with fasta_path.open("r") as fasta_reader:
            for i, contig in enumerate(SeqIO.parse(fasta_reader, "fasta")):
                # Skip contigs below Prodigal minimum length
                if len(contig) < MIN_LENGTH:
                    # print(f"contig {i + 1} from {str(fasta_path)} is too short ({len(contig)})")
                    continue

                predict_cds_features(contig, meta=len(contig) < META_LENGTH)
                dataframe = contig_to_dataframe(contig)
                dfs.append(dataframe)

    return dfs


def main(arguments):
    """
    Main function to interface between command line and training function.

    :param arguments: command line arguments
    :type arguments: list
    """
    args = parse_args(arguments)

    prophage_dir = args.prophages
    if not prophage_dir.is_dir():
        print(f"'{str(prophage_dir)}' is not a valid prophage directory")
        sys.exit(1)

    bacteria_dir = args.bacteria
    if not bacteria_dir.is_dir():
        print(f"'{str(bacteria_dir)}' is not a valid bacteria directory")
        sys.exit(1)

    print("Annotating prophage genes, and analyzing features...")
    prophage_dfs = annotate_sequences(prophage_dir)
    prophage_df = pd.concat(prophage_dfs, axis=0)
    prophage_df["category"] = [1] * len(prophage_df)
    # prophage_df.to_csv("prophages.csv", index=False)
    # prophage_df.corr().to_csv("prophage_correlation_matrix.csv")

    print("Annotating bacterial genes, and analyzing features...")
    bacteria_dfs = annotate_sequences(bacteria_dir)
    bacteria_df = pd.concat(bacteria_dfs, axis=0)
    bacteria_df["category"] = [0] * len(bacteria_df)
    # bacteria_df.to_csv("bacteria.csv", index=False)
    # bacteria_df.corr().to_csv("bacteria_correlation_matrix.csv")

    print("Training RandomForestClassifier...")
    new_model = train_random_forest_classifier(prophage_df, bacteria_df)
    print("Best model feature weights:")
    for feature_name, feature_weight in zip(prophage_df.columns[:-1],
                                            new_model.feature_importances_):
        print(f"{feature_name: <11}: {feature_weight * 100:.2f}%")

    # Save old model file if it exists
    if MODEL_PATH.is_file():
        old_path = MODEL_PATH.with_name("old.prophage_model.pickle")
        with MODEL_PATH.open("rb") as model_reader:
            old_model = pickle.load(model_reader)
        with old_path.open("wb") as model_writer:
            pickle.dump(old_model, model_writer)

    # Now save new model
    with MODEL_PATH.open("wb") as model_writer:
        pickle.dump(new_model, model_writer)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main(sys.argv[1:])
