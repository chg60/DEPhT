import argparse
import time

from prophicient_utils.pipelines import pull_sequences


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
PIPELINES = ["build_reference_db", "build_functions_db", "pull_sequences"]


def main(unparsed_args):
    args = parse_prophicient_utilities(unparsed_args)  

    start = time.time()

    if args.pipeline == "build_reference_db":
        pass
    elif args.pipeline == "build_functions_db":
        pass
    elif args.pipeline == "pull_sequences":
        pull_sequences.main(unparsed_args[1:])
    else:
        raise NotImplementedError(
                   f"Prophicient Utility pipeline '{args.build_reference_db}' "
                    "is not supported.")

    stop = time.time()

    print("\n\nPipeline completed.\nTime elapsed {:.2f}".format(stop - start))


def parse_prophicient_utilities(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("pipeline", type=str, choices=PIPELINES)

    args = parser.parse_args(unparsed_args[:1])

    return args
