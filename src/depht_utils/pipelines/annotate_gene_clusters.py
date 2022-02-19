import argparse
import pathlib
import sys

from Bio import SeqIO

from depht_utils.data.defaults import HHSUITEDB_DEFAULTS as DEFAULTS
from depht_utils.functions.fileio import (
    read_gene_index_file,
    write_cluster_function_index_file)

# GLOBAL VARIABLES
DEFAULTS = {"name": DEFAULTS["name"],
            "default_product": DEFAULTS["default_product"],
            "cutoff": DEFAULTS["consensus_annotation_cutoff"]}


# MAIN FUNCTIONS
def parse_annotate_gene_clusters(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("fasta_dir", type=pathlib.Path)
    parser.add_argument("index_file", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--name", type=str)

    args = parser.parse_args(unparsed_args)
    return args


def annotate_gene_clusters(fasta_dir, index_file, output_dir,
                           name=DEFAULTS["name"],
                           default_product=DEFAULTS["default_product"],
                           cutoff=DEFAULTS["cutoff"]):
    gene_index = read_gene_index_file(index_file)

    gene_cluster_data = []
    for input_file in fasta_dir.iterdir():
        records = [record for record in SeqIO.parse(input_file, "fasta")]

        product_counts = {}
        for record in records:
            index = record.id

            gene_data = gene_index.get(index)
            if gene_data is None:
                continue

            product_count = product_counts.get(gene_data["product"], 0)
            product_count += 1
            product_counts[gene_data["product"]] = product_count

        products = [product for product, count in product_counts.items()
                    if (count / len(records)) >= cutoff]
        products.sort(key=lambda x: product_counts[x], reverse=True)

        cluster_product = default_product
        for product in products:
            if product == default_product:
                continue

            cluster_product = product
            break

        gene_cluster_data.append([input_file.name, cluster_product])

    output_dir.mkdir(exist_ok=True, parents=True)
    out_index_file = output_dir.joinpath(".".join([name, "ppi"]))

    write_cluster_function_index_file(gene_cluster_data, out_index_file)

    return out_index_file


def main(unparsed_args):
    args = parse_annotate_gene_clusters(unparsed_args)
    annotate_gene_clusters(args.fasta_dir, args.index_file, args.output_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
