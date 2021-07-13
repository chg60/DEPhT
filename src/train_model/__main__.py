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
import matplotlib.pyplot as plt

from prophicient.functions.prophage_prediction import MODEL_PATH
from prophicient.functions.annotation import annotate_contig
from prophicient.functions.prophage_prediction import calculate_feature_dict
from train_model.functions.model_training import train_prophage_classifier

EPILOG = """
Bacterial genomes used for training should be complete assemblies, free of 
prophages, and with all plasmid contigs removed.
"""
TMP_DIR = pathlib.Path("/tmp/train_model")
# MODEL_PATH = pathlib.Path("~/Desktop/prophage_model.pickle").expanduser()


def parse_args(arguments):
    """
    Uses argparse module to parse the command line arguments.

    :param arguments: command line arguments
    :type arguments: list
    :return: parsed_args
    """
    p = argparse.ArgumentParser(description=__doc__, epilog=EPILOG)
    p.add_argument("-p", "--phage-dir", type=pathlib.Path, required=True,
                   help="path to a directory containing (pro)phage FASTA files")
    p.add_argument("-b", "--bacteria-dir", type=pathlib.Path, required=True,
                   help="path to a directory containing bacterial FASTA files")
    p.add_argument("--tmp-dir", type=pathlib.Path, default=TMP_DIR,
                   help=f"temporary directory for file I/O [default: "
                        f"{TMP_DIR}]")
    return p.parse_args(arguments)


def get_dataframe(filepath, tmp_dir):
    dataframes = list()

    for f in filepath.iterdir():
        if f.suffix not in (".fasta", ".fna"):
            continue
        record = SeqIO.read(f, "fasta")
        annotate_contig(record, tmp_dir, trna=False)

        starts, stops, strands, temp_dict = list(), list(), list(), dict()
        for feature in record.features:
            if feature.type == "CDS":
                starts.append(feature.location.start)
                stops.append(feature.location.end)
                strands.append(feature.location.strand)

        fd = calculate_feature_dict(record)

        temp_dict["contig_id"] = [record.id] * len(starts)
        temp_dict["start"] = starts
        temp_dict["stop"] = stops
        temp_dict["strand"] = strands
        temp_dict["lag_size"] = fd["lagging_window"]["gene_size"]
        temp_dict["ctr_size"] = fd["center_window"]["gene_size"]
        temp_dict["lead_size"] = fd["leading_window"]["gene_size"]
        temp_dict["lag_strand"] = fd["lagging_window"]["strand_change"]
        temp_dict["ctr_strand"] = fd["center_window"]["strand_change"]
        temp_dict["lead_strand"] = fd["leading_window"]["strand_change"]

        dataframes.append(pd.DataFrame(temp_dict))

    return pd.concat([pd.DataFrame(fd) for fd in dataframes], axis=0)


def main(arguments):
    """
    Train a model on the input data.

    :param arguments: command line arguments
    :type arguments: list
    """
    args = parse_args(arguments)

    # Make sure prophage directory exists
    prophage_dir = args.phage_dir
    if not prophage_dir.is_dir():
        print(f"'{str(prophage_dir)}' is not a valid prophage directory")
        sys.exit(1)

    # Make sure bacteria directory exists
    bacteria_dir = args.bacteria_dir
    if not bacteria_dir.is_dir():
        print(f"'{str(bacteria_dir)}' is not a valid bacteria directory")
        sys.exit(1)

    # Refresh temporary directory if it exists
    tmp_dir = args.tmp_dir
    if tmp_dir.is_dir():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)

    # Deal with the phage data
    prophage_file = prophage_dir.joinpath("data.csv")
    if prophage_file.is_file():
        print("Using existing prophage data...")
        prophage_df = pd.read_csv(prophage_file)
    else:
        print("Annotating and analyzing (pro)phages...")
        prophage_df = get_dataframe(prophage_dir, tmp_dir)
        prophage_df["is_prophage"] = [1] * len(prophage_df)
        prophage_df.to_csv(prophage_file, index=False)
    prophage_df = prophage_df.loc[:, ["ctr_size", "ctr_strand", "is_prophage"]]

    # prophage_df.plot(x="gene_size", y="strand_change", kind="scatter")
    # plt.show()

    # Deal with bacteria data
    bacteria_file = bacteria_dir.joinpath("data.csv")
    if bacteria_file.is_file():
        print("Using existing bacterial data...")
        bacteria_df = pd.read_csv(bacteria_file)
    else:
        print("Annotating and analyzing bacteria...")
        bacteria_df = get_dataframe(bacteria_dir, tmp_dir)
        bacteria_df["is_prophage"] = [0] * len(bacteria_df)
        bacteria_df.to_csv(bacteria_file, index=False)
    bacteria_df = bacteria_df.loc[:, ["ctr_size", "ctr_strand", "is_prophage"]]

    # bacteria_df.plot(x="gene_size", y="strand_change", kind="scatter")
    # plt.show()

    print("Training classifier...")
    model, training_mcc = train_prophage_classifier(prophage_df, bacteria_df)
    print(f"ProphageClassifier got average MCC = {training_mcc:.3f}")

    # Save old model file if it exists
    if MODEL_PATH.is_file():
        model_path = MODEL_PATH.with_name("new.prophage_model.pickle")
    else:
        model_path = MODEL_PATH

    # Now save model
    with model_path.open("wb") as model_writer:
        pickle.dump(model, model_writer)

    # Clean up before we finish
    shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main(sys.argv[1:])
