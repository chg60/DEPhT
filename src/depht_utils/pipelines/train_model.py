"""
Template for new model training pipeline.
"""

import sys
import argparse
import pathlib


def parse_args(unparsed_args=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-b", "--bact-dir", type=pathlib.Path,
                   help="path to directory containing bacterial file(s) for "
                        "training")
    p.add_argument("-p", "--phage-dir", type=pathlib.Path,
                   help="path to directory containing phage file(s) for "
                        "training")
    p.add_argument("-n", "--model-name", type=str,
                   help="name for the new model (e.g. 'Streptomyces')")

    if unparsed_args:
        return p.parse_args(unparsed_args)
    else:
        return p.parse_args()


def main(unparsed_args=None):
    args = parse_args(unparsed_args)
    print(args)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main()
