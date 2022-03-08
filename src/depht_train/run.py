"""Program to run DEPhT training pipelines."""

import argparse
import sys
import time

from depht_train.pipelines import (
    benchmark_output, build_reference_db, build_HMM_db, create_model,
    curate_gene_clusters, index_sequences, phamerate,
    pull_sequences, screen_conserved_phams, train_model)

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
PIPELINES = [
    "benchmark_output", "build_reference_db", "build_HMM_db", "create_model",
    "curate_gene_clusters", "index_sequences", "phamerate", "pull_sequences",
    "screen_conserved_phams", "train_model"]

EPILOG = "The 'create_model' pipeline performs all model-building steps. " \
         "The other pipelines can be used in combination for better control " \
         "over the resultant model."


def parse_args(unparsed_args):
    """Use argparse to verify pipeline argument only.

    :param unparsed_args:
    :type unparsed_args:
    """
    parser = argparse.ArgumentParser(description=__doc__, prog="depht_train",
                                     epilog=EPILOG)

    parser.add_argument("pipeline", type=str, choices=PIPELINES,
                        metavar="pipeline",
                        help=f"name of the depht_train pipeline to run "
                             f"[choices: {PIPELINES}]")

    return parser.parse_args(unparsed_args[1:2])


def main():
    """Commandline entrypoint to this module."""
    if len(sys.argv) == 1:
        sys.argv.append("-h")

    unparsed_args = sys.argv

    args = parse_args(unparsed_args)

    start = time.time()

    if args.pipeline == "benchmark_output":
        benchmark_output.main(unparsed_args[1:])
    elif args.pipeline == "build_reference_db":
        build_reference_db.main(unparsed_args[1:])
    elif args.pipeline == "build_HMM_db":
        build_HMM_db.main(unparsed_args[1:])
    elif args.pipeline == "create_model":
        create_model.main(unparsed_args[1:])
    elif args.pipeline == "curate_gene_clusters":
        curate_gene_clusters.main(unparsed_args[1:])
    elif args.pipeline == "index_sequences":
        index_sequences.main(unparsed_args[1:])
    elif args.pipeline == "phamerate":
        phamerate.main(unparsed_args[1:])
    elif args.pipeline == "pull_sequences":
        pull_sequences.main(unparsed_args[1:])
    elif args.pipeline == "screen_conserved_phams":
        screen_conserved_phams.main(unparsed_args[1:])
    elif args.pipeline == "train_model":
        train_model.main(unparsed_args[1:])
    else:
        raise NotImplementedError(f"unknown depht_train pipeline "
                                  f"{args.pipeline}")

    stop = time.time()

    print("\n\nPipeline completed.\nTime elapsed {:.2f}s".format(stop - start))
