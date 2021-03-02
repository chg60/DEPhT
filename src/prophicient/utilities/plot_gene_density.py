import sys
import argparse
import pathlib
import matplotlib.pyplot as plt

from Bio.SeqIO import read

BIN_WIDTH = 1000
WINDOW_SIZE = 50000


def parse_args(arguments):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("genbank_directory", type=pathlib.Path, help="directory containing Genbank flatfiles")
    return p.parse_args(arguments)


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

    # Initialize every interval's counter to 0
    for interval_index in range(0, len(record), interval):
        genes_per_interval[interval_index] = 0

    # Iterate over genes, and place them in their appropriate bin
    for gene in record.features:
        interval_index = gene.location.start // interval * interval
        genes_per_interval[interval_index] += 1

    return genes_per_interval


def plot_gene_density(interval_counts, record_length, average_density, name, window_size=WINDOW_SIZE, advance=BIN_WIDTH):
    # Now calculate xs and ys for scatterplot
    xs, ys, zs = list(), list(), list()
    for x in range(0, record_length - window_size, advance):
        total = 0
        for j in range(x, x + window_size, advance):
            total += interval_counts[j]
        y = total / window_size * advance
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
    plt.show()


def main(arguments):
    args = parse_args(arguments)
    gen_dir = args.genbank_directory

    # Make sure Genbank directory exists...
    if not gen_dir.is_dir():
        print(f"'{str(gen_dir)}' is not a valid directory")
        sys.exit(1)

    # Walk Genbank directory to get flatfile names
    for flatfile_path in gen_dir.iterdir():
        with flatfile_path.open("r") as fh:
            try:
                record = read(fh, "genbank")
            except ValueError:
                continue

            avg_dens = float(len(record.features)) / len(record) * BIN_WIDTH
            interval_genes = count_genes_per_interval(record)
            plot_gene_density(interval_genes, len(record), avg_dens, flatfile_path.stem)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main(sys.argv[1:])
