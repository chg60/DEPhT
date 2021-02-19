"""Pipeline for downloading, clustering, and making hhsuite3-compatible
databases from phage genomes."""
import argparse
import configparser
from configparser import Error as ConfigError
from pathlib import Path
import sys

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_CONFIG_SECTIONS = {"ncbi": {"api_key", "email", "tool"}}
DEFAULT_FOLDER_NAME = "PhageGeneRef"


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the compile db pipeline

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    args = parse_compile_db(unparsed_args_list)

    config = build_complete_config(args.config_file)
    config

    taxon_ids = parse_taxon_ids(args.identifiers_file)

    execute_compile_db(taxon_ids, folder_path=args.folder_path,
                       folder_name=args.folder_name, verbose=args.verbose,
                       force=args.force, cores=args.number_processes)
    pass


def parse_compile_db(unparsed_args_list):
    """Parses export_db arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :typ unparsed_args_list: list[str]
    :returns: Argparse module parsed args.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("identifiers_file", type=convert_file_path)

    parser.add_argument("-m", "--folder_name", type=str)
    parser.add_argument("-o", "--folder_path", type=Path)
    parser.add_argument("-c", "--config_file", type=convert_file_path)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-f", "--force", action="store_true")
    parser.add_argument("-np", "--number_processes", type=int)

    parser.set_defaults(folder_path=None, folder_name=DEFAULT_FOLDER_NAME,
                        config_file=None, verbose=False, force=False,
                        number_processes=1)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_compile_db(taxon_ids, folder_path=None,
                       folder_name=DEFAULT_FOLDER_NAME, verbose=False,
                       force=False, cores=1):
    pass


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def convert_file_path(path_string):
    """Function to convert argparse input to a working file path.

    :param path: A string to be converted into a Path object.
    :type path: str
    :returns: A Path object converted from the inputted string.
    :rtype: pathlib.Path
    """
    path = Path(path_string)
    path = path.expanduser()
    path = path.resolve()

    if not path.is_file():
        print(f"The file {path} does not exist.")
        sys.exit(1)

    return path


def parse_taxon_ids(identifiers_file):
    """Parses a taxon identifiers file, assuming each taxon ID is on
    a separate line.

    :param identifiers_file: Path to a file containing taxon identifiers.
    :type identifiers_file: pathlib.Path
    :returns: A set of taxon identifiers.
    :rtype: set()
    """
    taxon_ids = set()

    with identifiers_file.open(mode="r") as filehandle:
        for line in filehandle:
            taxon_ids.add(line.rstrip())

    return taxon_ids


# CONFIG FILE SETUP
# -----------------------------------------------------------------------------
def setup_section(keys, value):
    dict = {}
    for key in keys:
        dict[key] = value
    return dict


def default_parser(null_value):
    """Constructs complete config with empty values."""
    # Need to allow no value if the null value is None.
    null_parser = configparser.ConfigParser(allow_no_value=True)
    config_sections = DEFAULT_CONFIG_SECTIONS

    for header in config_sections.keys():
        config_section = config_sections[header]
        for subheader in config_section.keys():
            config_section[subheader] = null_value

    return null_parser


def parse_config(config_path, parser=None):
    """Get parameters from config file."""
    if parser is None:
        parser = configparser.ConfigParser()

    try:
        parser.read(config_path)
    except ConfigError:
        print("Unable to parse config file.")
        sys.exit(1)
    else:
        return parser


def build_complete_config(config_path):
    "Buid a complete config object by merging user-supplied and default config"
    parser = default_parser(None)

    if config_path is not None:
        parser = parse_config(config_path, parser)

    return parser
