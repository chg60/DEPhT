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
import shutil

import pandas as pd
from Bio import SeqIO

from prophicient.classes.prophage import ANNOTATIONS
from prophicient.functions.annotation import annotate_contig
from prophicient.functions.prophage_prediction import calculate_feature_dict
from train_model.functions.model_training import train_random_forest_classifier

EPILOG = """
Bacterial genomes used for training should be complete assemblies, free of 
prophages, and with all plasmid contigs removed.
"""
TMP_DIR = pathlib.Path("/tmp/train_model")
MODEL_PATH = pathlib.Path("~/Desktop/prophage_model.pickle").expanduser()


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
    p.add_argument("--tmp-dir", type=pathlib.Path, default=TMP_DIR,
                   help=f"temporary directory for file I/O [default: "
                        f"{TMP_DIR}]")
    return p.parse_args(arguments)


def fasta_to_gbk(fasta_file, tmp_dir):
    """
    Annotates CDS features on the sequence in the indicated fasta file.

    Returns the path to an annotated Genbank record for this sequence.

    :param fasta_file: the path to a FASTA nucleotide sequence file
    :type fasta_file: pathlib.Path
    :param tmp_dir: path where temporary files can go
    :type tmp_dir: pathlib.Path
    :return: genbank_file
    """
    genbank_file = fasta_file.with_suffix(".gbk")

    if not genbank_file.is_file():
        record = SeqIO.read(fasta_file, "fasta")
        annotate_contig(record, tmp_dir)
        record.annotations = ANNOTATIONS
        SeqIO.write(record, genbank_file, "genbank")
    return genbank_file


def gbk_to_csv(genbank_file):
    """
    Builds the feature dataframe from the SeqRecord in a Genbank
    flatfile. Dumps the dataframe to a CSV file, and returns that
    filepath

    :param genbank_file: the path to a Genbank flat file
    :type genbank_file: pathlib.Path
    :return: dataframe_file
    """
    dataframe_file = genbank_file.with_suffix(".csv")

    if not dataframe_file.is_file():
        record = SeqIO.read(genbank_file, "genbank")
        df = pd.DataFrame(calculate_feature_dict(record))
        df.to_csv(dataframe_file, index=False)
    return dataframe_file


def get_dataframe(filepath, tmp_dir):
    dataframes = list()
    for f in filepath.iterdir():
        if f.suffix not in (".fasta", ".fna"):
            continue
        record = SeqIO.read(f, "fasta")
        annotate_contig(record, tmp_dir, no_trna=True)
        dataframes.append(pd.DataFrame(calculate_feature_dict(record)))
    return pd.concat([pd.DataFrame(fd) for fd in dataframes], axis=0)


def main(arguments):
    """
    Main function to interface between command line and training function.

    :param arguments: command line arguments
    :type arguments: list
    """
    args = parse_args(arguments)

    # Make sure prophage directory exists
    prophage_dir = args.prophages
    if not prophage_dir.is_dir():
        print(f"'{str(prophage_dir)}' is not a valid prophage directory")
        sys.exit(1)

    # Make sure bacteria directory exists
    bacteria_dir = args.bacteria
    if not bacteria_dir.is_dir():
        print(f"'{str(bacteria_dir)}' is not a valid bacteria directory")
        sys.exit(1)

    # Make sure temporary directory exists
    tmp_dir = args.tmp_dir
    if not tmp_dir.is_dir():
        tmp_dir.mkdir(parents=True)

    print("Annotating and analyzing prophages...")
    prophage_df = get_dataframe(prophage_dir, tmp_dir)
    prophage_df["class"] = [1] * len(prophage_df)
    # prophage_df.to_csv("prophages.csv", index=False)
    # prophage_df.corr().to_csv("prophage_correlation_matrix.csv")

    print("Annotating and analyzing bacteria...")
    bacteria_df = get_dataframe(bacteria_dir, tmp_dir)
    bacteria_df["class"] = [0] * len(bacteria_df)
    # bacteria_df.to_csv("bacteria.csv", index=False)
    # bacteria_df.corr().to_csv("bacteria_correlation_matrix.csv")

    print("Training RandomForestClassifier...")
    new_model = train_random_forest_classifier(prophage_df, bacteria_df)
    print("Model feature weights:")
    for feature_name, feature_weight in zip(prophage_df.columns[:-1],
                                            new_model.feature_importances_):
        print(f"{feature_name: <11}: {feature_weight * 100:.2f}%")

    # Save old model file if it exists
    if MODEL_PATH.is_file():
        model_path = MODEL_PATH.with_name("new.prophage_model.pickle")
    else:
        model_path = MODEL_PATH

    # Now save new model
    with model_path.open("wb") as model_writer:
        pickle.dump(new_model, model_writer)

    # Clean up before we finish
    shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main(sys.argv[1:])
