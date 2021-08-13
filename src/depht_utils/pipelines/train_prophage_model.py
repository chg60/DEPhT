"""
Helper pipeline for training the model used by MycoPhinder on known
phage/prophage and bacterial gene data.
"""

import argparse
import pathlib
import pickle
import shutil
import sys

import pandas as pd
from Bio import SeqIO

from depht.classes.prophage_classifier import ProphageClassifier
from depht.functions.annotation import annotate_record
from depht.functions.prophage_prediction import MODEL_PATH, WINDOW
from depht.functions.prophage_prediction import build_contig_dataframe
from depht.functions.statistics import average, mcc, f1_score
from depht_utils.classes.kfold import KFold

EPILOG = """
Bacterial genomes used for training should be complete assemblies, free of 
prophages, and with all plasmid contigs removed.
"""
TMP_DIR = pathlib.Path("/tmp/prophicient_utils/train_model")


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
    p.add_argument("-t", "--tmp-dir", type=pathlib.Path, default=TMP_DIR,
                   help=f"temporary directory for file I/O [default: "
                        f"{TMP_DIR}]")
    p.add_argument("-w", "--window", type=int, default=WINDOW,
                   help=f"number of CDS genes to consider at once [default: "
                        f"{WINDOW}]")
    return p.parse_args(arguments)


def get_dataframe(filepath, tmp_dir, window):
    """
    Builds a bacterial/prophage dataframe from a directory of FASTA
    files.

    :param filepath: path to the FASTAs to integrate into a dataframe
    :type filepath: pathlib.Path
    :param tmp_dir: path where annotation files can be written
    :type tmp_dir: pathlib.Path
    :param window: how many genes to consider in a window
    :type window: int
    :return: dataframe
    """
    dataframes = list()

    for f in filepath.iterdir():
        if f.suffix not in (".fasta", ".fna"):
            continue
        record = SeqIO.read(f, "fasta")
        annotate_record(record, tmp_dir, trna=False)
        df = build_contig_dataframe(record, window)
        dataframes.append(df)

    return pd.concat(dataframes, axis=0)


def _score(real_classes, predict_classes):
    """
    Helper function to score class predictions and then return
    the mcc and f1_scores.
    """
    tps, tns, fps, fns = 0, 0, 0, 0

    for real_value, predict_value in zip(real_classes, predict_classes):
        if real_value and predict_value:
            tps += 1
        elif real_value and not predict_value:
            fns += 1
        elif predict_value and not real_value:
            fps += 1
        else:
            tns += 1

    print(tps, tns, fps, fns)

    return mcc(tps, fns, tns, fps), f1_score(tps, fps, fns)


def train_prophage_classifier(prophage_data, bacteria_data):
    """
    Trains a ProphageClassifier on the given prophage_data and
    bacteria_data. This is effectively a simple implementation of
    a Naive Bayes Classifier, that builds a probability distribution
    for each input feature, and can then use that distribution to
    predict the probability that a given input gene belongs to the
    prophage class.

    :param prophage_data: the (pro)phage training data
    :type prophage_data: pd.DataFrame
    :param bacteria_data: the bacterial training data
    :type bacteria_data: pd.DataFrame
    :return: clf, mean_mcc, mean_f1
    """
    # Seed RNG for reproducibility - use the "Answer to the Ultimate
    # Question of Life, the Universe, and Everything." (Douglas Adams)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    # Get train/test split indices for both groups of genes
    prophage_splits = [x for x in kf.split(len(prophage_data))]
    bacteria_splits = [x for x in kf.split(len(bacteria_data))]

    mcc_scores, f1_scores = list(), list()
    clf = ProphageClassifier()

    for prophage_split, bacteria_split in zip(prophage_splits, bacteria_splits):
        p_train, p_test = prophage_split
        b_train, b_test = bacteria_split

        train_feats = pd.concat((prophage_data.iloc[list(p_train), :-1],
                                 bacteria_data.iloc[list(b_train), :-1]),
                                axis=0)
        train_labels = pd.concat((prophage_data.iloc[list(p_train), -1],
                                  bacteria_data.iloc[list(b_train), -1]),
                                 axis=0)
        test_feats = pd.concat((prophage_data.iloc[list(p_test), :-1],
                                bacteria_data.iloc[list(b_test), :-1]), axis=0)
        test_labels = pd.concat((prophage_data.iloc[list(p_test), -1],
                                 bacteria_data.iloc[list(b_test), -1]), axis=0)

        clf.fit(train_feats, train_labels)

        predictions = clf.predict(test_feats, alpha=0.55)
        m, f = _score(test_labels, predictions)
        mcc_scores.append(m), f1_scores.append(f)

    mean_mcc, mean_f1 = average(mcc_scores), average(f1_scores)

    all_feats = pd.concat((prophage_data.iloc[:, :-1],
                           bacteria_data.iloc[:, :-1]), axis=0)
    all_labels = pd.concat((prophage_data.iloc[:, -1],
                            bacteria_data.iloc[:, -1]), axis=0)

    clf.fit(all_feats, all_labels, plot=True)

    return clf, mean_mcc, mean_f1


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

    # window size
    w = args.window

    # Deal with the phage data
    prophage_file = prophage_dir.joinpath(f"_{w}_data.csv")
    if prophage_file.is_file():
        print("Using existing prophage data...")
        prophage_df = pd.read_csv(prophage_file)
    else:
        print("Annotating and analyzing (pro)phages...")
        prophage_df = get_dataframe(prophage_dir, tmp_dir, w)
        prophage_df["is_prophage"] = [1] * len(prophage_df)
        prophage_df.to_csv(prophage_file, index=False)
    prophage_df = prophage_df.loc[:, ["ctr_size", "ctr_strand", "is_prophage"]]

    # Deal with bacteria data
    bacteria_file = bacteria_dir.joinpath(f"_{w}_data.csv")
    if bacteria_file.is_file():
        print("Using existing bacterial data...")
        bacteria_df = pd.read_csv(bacteria_file)
    else:
        print("Annotating and analyzing bacteria...")
        bacteria_df = get_dataframe(bacteria_dir, tmp_dir, w)
        bacteria_df["is_prophage"] = [0] * len(bacteria_df)
        bacteria_df.to_csv(bacteria_file, index=False)
    bacteria_df = bacteria_df.loc[:, ["ctr_size", "ctr_strand", "is_prophage"]]

    print("Training classifier...")
    model, train_mcc, train_f1 = \
        train_prophage_classifier(prophage_df, bacteria_df)
    print(f"ProphageClassifier got average MCC, F1 = {train_mcc:.3f}, "
          f"{train_f1:.3f}")

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
