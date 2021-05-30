import argparse
import pathlib
import sys

from Bio import SeqIO

from prophicient_utils.data.defaults import HHSUITEDB_DEFAULTS


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULTS = {"name": HHSUITEDB_DEFAULTS["name"], "table": 11,
            "default_product": HHSUITEDB_DEFAULTS["default_product"]}


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def parse_index_functions(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--name", type=str)

    parser.set_defaults(**DEFAULTS)
    args = parser.parse_args(unparsed_args)
    return args


def index_functions(input_dir, output_dir,
                    name=DEFAULTS["name"], table=DEFAULTS["table"],
                    default_product=DEFAULTS["default_product"]):
    cds_features = []
    for input_file in input_dir.iterdir():
        records = [record for record in SeqIO.parse(input_file, "gb")]

        for record in records:
            for index, feature in enumerate(record.features):
                if feature.type != "CDS":
                    continue

                # Checks to see if the CDS feature has an annotated translation
                translation = feature.qualifiers.get("translation")
                # If the feature lacks a translation qualifier, create one
                if not translation:
                    sequence = feature.extract(record.seq)
                    translation = sequence.translate(table=table, to_stop=True)

                    feature.qualifiers["translation"] = [translation]

                # Checks to see if the CDS feature has an annotated function
                product = feature.qualifiers.get("product")
                # If the feature lacks a function qualifier, create one
                if not product:
                    feature.qualifiers["product"] = [default_product]
                else:
                    # If the feature has a blank function qualifier, create one
                    if not product[0]:
                        feature.qualifiers["product"] = [default_product]

                # Checks to see if the CDS feature has an annotated locus tag
                locus_tag = feature.qualifiers.get("locus_tag")
                if not locus_tag:
                    feature.qualifiers["locus_tag"] = ["_".join(
                                                    [record.id, str(index+1)])]

                cds_features.append(feature)

    output_dir.mkdir(exist_ok=True, parents=True)

    fasta_file = output_dir.joinpath(".".join([name, "fasta"]))
    index_file = output_dir.joinpath(".".join([name, "pgi"]))

    with fasta_file.open(mode="w") as fasta_filehandle:
        with index_file.open(mode="w") as index_filehandle:
            for index, feature in enumerate(cds_features):
                locus_tag = feature.qualifiers["locus_tag"][0]
                translation = feature.qualifiers["translation"][0]
                product = feature.qualifiers["product"][0]

                fasta_filehandle.write("".join([">", str(index), "\n"]))
                fasta_filehandle.write("".join([translation, "\n"]))

                index_filehandle.write("\t".join(
                                       [str(index), locus_tag, product, "\n"]))


def main(unparsed_args):
    args = parse_index_functions(unparsed_args)
    index_functions(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
