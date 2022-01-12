import argparse
import time

from depht_utils.pipelines import (
    benchmark_output, build_reference_db,
    build_HMM_db, curate_gene_clusters, index_sequences, phamerate,
    pull_sequences, screen_conserved_phams)

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
PIPELINES = [
        "benchmark_output", "build_reference_db",
        "build_HMM_db", "curate_gene_clusters", "index_sequences",
        "phamerate", "pull_sequences", "screen_conserved_phams",
        "train_prophage_model", "train_model"]


def main(unparsed_args):
    args = parse_prophicient_utilities(unparsed_args)

    start = time.time()

    if args.pipeline == "benchmark_output":
        benchmark_output.main(unparsed_args[1:])
    elif args.pipeline == "build_reference_db":
        build_reference_db.main(unparsed_args[1:])
    elif args.pipeline == "build_HMM_db":
        build_HMM_db.main(unparsed_args[1:])
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
    elif args.pipeline == "train_prophage_model":
        train_prophage_model.main(unparsed_args[1:])
    elif args.pipeline == "train_model":
        train_model.main(unparsed_args[1:])
    else:
        raise NotImplementedError(
                   f"Prophicient Utility pipeline '{args.build_reference_db}' "
                   "is not supported.")

    stop = time.time()

    print("\n\nPipeline completed.\nTime elapsed {:.2f}s".format(stop - start))


def parse_prophicient_utilities(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("pipeline", type=str, choices=PIPELINES)

    args = parser.parse_args(unparsed_args[:1])

    return args
