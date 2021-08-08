import configparser
import sys
from configparser import Error as ConfigError

DEFAULT_CONFIG_SECTIONS = {"ncbi": {"api_key", "email", "tool"}}


def setup_section(keys, value):
    section = {}
    for key in keys:
        section[key] = value
    return section


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
