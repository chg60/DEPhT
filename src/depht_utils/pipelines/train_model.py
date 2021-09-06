"""
Template for new model training pipeline.
"""

import sys
import argparse
import pathlib


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    return p.parse_args()


def main():
    pass


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main()
