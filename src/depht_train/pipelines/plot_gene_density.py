"""Utility script to plot gene density along a genome.
"""
import argparse
import pathlib
import sys

import matplotlib.pyplot as plt
from Bio import SeqIO

# Defaults for bin width and window size, in kb
BIN_WIDTH = 1000
WINDOW_SIZE = 50000


def parse_args(unparsed_args):
    """Parse commandline arguments.

    :param unparsed_args:
    :type unparsed_args: list[str]
    """
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("flatfile_dir", type=pathlib.Path,
                   help="path to a directory of Genbank flatfiles to inspect")

    return p.parse_args(unparsed_args)


def count_genes_per_interval(record, interval=BIN_WIDTH):
    """
    Function to count the number of genes per `interval` nucleotides in the record.

    Counts gene starts, but could easily be adapted to work in some other way

    :param record: a Genbank record to process
    :type record: Bio.SeqRecord.SeqRecord
    :param interval: the number of bases to advance when indexing
    :type interval: int
    :return: genes_per_interval
    """

    genes_per_interval = dict()
    max_interval = len(record) // interval * interval

    # Initialize every interval's counter to 0
    for interval_index in range(0, len(record), interval):
        genes_per_interval[interval_index] = 0

    # Iterate over genes, and place them in their appropriate bin
    for gene in record.features:
        if gene.type == "CDS":
            interval_index = gene.location.start // interval * interval
            genes_per_interval[interval_index] += 1

    # Wrap around ends of the genome
    for interval_index in range(0, WINDOW_SIZE, interval):
        genes_per_interval[interval_index + max_interval] = genes_per_interval[interval_index]
    for interval_index in range(max_interval - WINDOW_SIZE, max_interval, interval):
        genes_per_interval[interval_index - max_interval] = genes_per_interval[interval_index]

    return genes_per_interval


def create_gene_density_plot(interval_counts, record_length, average_density, name, window_size=WINDOW_SIZE, advance=BIN_WIDTH):
    # Now calculate xs and ys for scatterplot
    xs, ys, zs = list(), list(), list()
    for x in range(0, record_length, advance):
        total = 0
        for j in range(x - window_size // 2, x + window_size // 2, advance):
            total += interval_counts[j]
        y = (total / window_size) * advance
        z = (y - average_density) ** 2
        xs.append(x)
        ys.append(y)
        zs.append(z)
    stdev = (float(sum(zs))/len(zs)) ** 0.5
    plt.scatter(xs, ys)
    plt.axhline(y=average_density, color="r", linestyle="-")
    plt.axhline(y=average_density + 2*stdev, color="r", linestyle="--")
    plt.axhline(y=average_density - 2*stdev, color="r", linestyle="--")
    plt.title(f"{name} Gene Density")
    plt.xlabel("Genomic Position (Mbp)")
    plt.ylabel(f"genes/kb over {window_size}kb window")
    plt.ylim(0.25, 2.25)
    return max(ys)


def get_records(flatfile):
    """
    Parses and returns the record from the indicated Genbank flatfile.

    :param flatfile: the path to the flatfile to fetch record(s) from
    :type flatfile: pathlib.Path
    :return: records
    """
    records = list()
    with flatfile.open("r") as fh:
        for record in SeqIO.parse(fh, "genbank"):
            records.append(record)
    return records


def main(unparsed_args=None):
    """Commandline entrypoint to this module."""
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])
    ff_dir = args.flatfile_dir

    # Make sure Genbank directory exists...
    if not ff_dir.is_dir():
        print(f"input directory '{str(ff_dir)}' does not exist... exiting")
        sys.exit(1)

    # Iterate over ff_dir
    for flatfile in ff_dir.iterdir():
        print(f"Processing {str(flatfile)}...")
        # Parse the Genbank record
        try:
            records = get_records(flatfile)
        except:
            print(f"'{str(flatfile)} is fucked.")
            continue

        if len(records) > 5:
            print(f"'{str(flatfile)}' contains {len(records)} records... is this assembly complete?")

        for record in records:
            if len(record) < WINDOW_SIZE:
                print(f"\tcontig '{record.id}' is too short to inspect for prophages...")
                continue

            avg_dens = float(len([feature for feature in record.features if feature.type == "CDS"])) / len(record) * BIN_WIDTH
            interval_genes = count_genes_per_interval(record)
            max_value = create_gene_density_plot(interval_genes, len(record), avg_dens, flatfile.stem)
            if max_value >= 1.5:
                save_path = flatfile.with_name(f"{flatfile.stem}_x.png")
            else:
                save_path = flatfile.with_suffix(".png")
            plt.savefig(save_path)
            plt.close()


if __name__ == "__main__":
    main()
