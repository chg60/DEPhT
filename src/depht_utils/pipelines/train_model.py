"""Standalone script for training just the phage/bacterial classifier
portion of the model."""

import sys
import argparse
import pathlib
import pickle

import pandas as pd
from Bio import SeqIO

from depht.functions.annotation import (annotate_record,
                                        cleanup_flatfile_records)
from depht.functions.multiprocess import parallelize, LOGICAL_CORES
from depht.functions.prophage_prediction import build_contig_dataframe
from depht.functions.sniff_format import sniff_format
from depht_utils.functions.train_classifier import train_classifier


MODEL_DIR = pathlib.Path().home().joinpath(".depht/models")
WINDOW = 55


def parse_args(unparsed_args):
    """Parse commandline arguments."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("name", type=str,
                   help="one-word name for the resultant model")
    p.add_argument("phages", type=pathlib.Path,
                   help="path to a directory containing phage sequences")
    p.add_argument("bacteria", type=pathlib.Path,
                   help="path to a directory containing bacterial sequences")
    p.add_argument("-p", "--prophage-csv", type=pathlib.Path, default=None,
                   help="path to a 3-column CSV mapping the coordinates of "
                        "any known prophages in the bacterial contigs")
    p.add_argument("-w", "--window-size", type=int, default=WINDOW,
                   help=f"number of genes to average features over "
                        f"[default: {WINDOW}]")
    p.add_argument("-c", "--cpu-cores", type=int, default=LOGICAL_CORES,
                   help=f"number of cpu cores to use "
                        f"[default: {LOGICAL_CORES}]")
    return p.parse_args(unparsed_args)


def get_genome_df(filepath, tmp_dir, window):
    """
    Builds a single bacterial/prophage genome's dataframe from a
    directory of FASTA files.

    :param filepath:
    :type filepath:
    :param tmp_dir:
    :type tmp_dir:
    :param window:
    :type window:
    :return: df
    """
    df_list = list()

    file_fmt = sniff_format(filepath)
    records = sorted([x for x in SeqIO.parse(filepath, file_fmt)],
                     reverse=True, key=lambda x: len(x))

    if file_fmt == "fasta":
        for record in records:
            annotate_record(record, tmp_dir, trna=False)
    else:
        cleanup_flatfile_records(records)

    for record in records:
        if len(record.features) > window:
            df_list.append(build_contig_dataframe(record, window))

    if len(df_list) > 0:
        return pd.concat(df_list, axis=0)


def get_dataframe(filepath, tmp_dir, window, cpus=1):
    """
    Builds a complete bacterial/prophage dataframe from a directory
    of FASTA files.

    :param filepath: path to the FASTAs to integrate into a dataframe
    :type filepath: pathlib.Path
    :param tmp_dir: path where annotation files can be written
    :type tmp_dir: pathlib.Path
    :param window: how many genes to consider in a window
    :type window: int
    :param cpus:
    :type cpus:
    :return: dataframe
    """
    if not tmp_dir.is_dir():
        tmp_dir.mkdir(parents=True)

    jobs = list()
    for f in filepath.iterdir():
        file_fmt = sniff_format(f)

        if file_fmt not in ("fasta", "genbank"):
            continue

        jobs.append((f, tmp_dir, window))

    dataframes = parallelize(jobs, cpus, get_genome_df)

    return pd.concat(dataframes, axis=0)


def main(unparsed_args):
    """Commandline entry point for the script."""
    namespace = parse_args(unparsed_args)

    name = namespace.name
    phg_dir = namespace.phages
    bct_dir = namespace.bacteria
    window = namespace.window_size
    prophages = namespace.prophage_csv
    cpus = namespace.cpu_cores

    train_model(name, phg_dir, bct_dir, window=window, prophages=prophages,
                cpus=cpus)


def train_model(name, phg_dir, bct_dir,
                window=WINDOW, prophages=None, cpus=LOGICAL_CORES):
    if not phg_dir.is_dir():
        print(f"specified phage directory '{str(phg_dir)}' is not a valid "
              f"directory")
        return

    if not bct_dir.is_dir():
        print(f"specified bacteria directory '{str(bct_dir)}' is not a valid "
              f"directory")
        return

    model_dir = MODEL_DIR.joinpath(name)

    tmp_dir = model_dir.joinpath("tmp")
    tmp_dir.mkdir(exist_ok=True, parents=True)

    phg_csv = model_dir.joinpath(f"_{window}_phage_df.csv")
    if phg_csv.is_file():
        phg_df = pd.read_csv(phg_csv)
    else:
        phg_df = get_dataframe(phg_dir, tmp_dir, window, cpus)
        phg_df["class"] = [1] * len(phg_df)
        phg_df.to_csv(phg_csv, index=False)

    bct_csv = model_dir.joinpath(f"_{window}_bacteria_df.csv")
    if bct_csv.is_file():
        bct_df = pd.read_csv(bct_csv)
    else:
        bct_df = get_dataframe(bct_dir, tmp_dir, window, cpus)
        bct_df["class"] = [0] * len(bct_df)

        if prophages:
            # prevent SettingWithCopyWarning message from appearing
            pd.options.mode.chained_assignment = None

            pro_df = pd.read_csv(prophages, dtype={"contig_id": str,
                                                   "start": int,
                                                   "stop": int})
            for i, prophage in pro_df.iterrows():
                name = prophage["contig_id"]
                start, stop = prophage["start"], prophage["stop"]

                bct_df["class"].loc[
                    (bct_df["contig_id"] == name) &
                    (bct_df["start"] >= start) &
                    (bct_df["stop"] >= start) &
                    (bct_df["start"] <= stop) &
                    (bct_df["stop"] <= stop)] = 1

        bct_df.to_csv(bct_csv, index=False)

    phg_df = phg_df.loc[:, ["ctr_size", "ctr_strand", "class"]]
    bct_df = bct_df.loc[:, ["ctr_size", "ctr_strand", "class"]]

    joint_df = pd.concat([bct_df, phg_df], axis=0)

    clf, figs = train_classifier(joint_df, verbose=True)
    for feature_name, fig in figs:
        filename = model_dir.joinpath(f"{name}_{feature_name}_w{window}.pdf")
        fig.update_layout(title=f"{name}: {feature_name} (w={window})")
        if feature_name == "ctr_size":
            fig.update_xaxes(range=[0, 2000])
        else:
            fig.update_xaxes(range=[0, 35])
        fig.write_image(filename)

    clf_path = model_dir.joinpath("classifier.pkl")
    with open(clf_path, "wb") as classifier_writer:
        pickle.dump(clf, classifier_writer)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main()
