"""Helper script to facilitate downloading models from OSF.io."""

import argparse
import pathlib
import requests
import shutil
import sys
import zipfile

DEPHT_DIR = pathlib.Path().home().joinpath(".depht/models")

URLS = {"Gordonia": "https://osf.io/download/djwsb/",
        "Mycobacterium": "https://osf.io/download/aw4up/"}


def parse_args():
    """Parse command line arguments."""
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("model",
                   type=str, choices=["Mycobacterium", "Gordonia"],
                   help="choose a model to download")

    p.add_argument("-f", "--force", action="store_true",
                   help="overwrite existing model files if they exist "
                        "[default: False]")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="print verbose output to stdout [default: False]")

    return p.parse_args()


def download_model(url, path):
    """Download model from OSF.io and save to path."""
    r = requests.get(url, allow_redirects=True)
    open(path, "wb").write(r.content)

    return path


def unzip_model(path):
    """Unzip model to path."""
    with zipfile.ZipFile(path, "r") as zip_ref:
        zip_ref.extractall(path.parent)


def main():
    """Commandline entry point."""
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = parse_args()

    model = args.model
    model_url = URLS[model]
    zip_path = DEPHT_DIR.joinpath(f"{model}.zip")
    model_path = DEPHT_DIR.joinpath(model)
    overwrite = args.force

    # Create the DEPhT data directory if it doesn't exist
    if not DEPHT_DIR.is_dir():
        DEPHT_DIR.mkdir(parents=True)

    # If model exists and not told to overwrite, exit
    if model_path.is_dir() and not overwrite:
        print(f"Model files already exist at {model_path}. "
              "Use '-f/--force' to overwrite.")
        sys.exit(1)

    # If model exists and told to overwrite, remove it
    if model_path.is_dir() and overwrite:
        if args.verbose:
            print(f"Removing existing model files at {model_path}...")
        shutil.rmtree(model_path)

    # Download and unzip model
    if args.verbose:
        print(f"Downloading model files from {model_url}...")
    download_model(model_url, zip_path)
    if args.verbose:
        print(f"Unzipping model files to {model_path}...")
    unzip_model(zip_path)
    if args.verbose:
        print(f"Removing zip file at {zip_path}...")
    zip_path.unlink()
    if args.verbose:
        print("Done.")


if __name__ == "__main__":
    main()
