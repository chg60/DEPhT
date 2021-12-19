"""Function that tries to detect biological file format by sniffing the
first line of a given file."""


def sniff_format(filepath):
    """Read the first line of the file and return its predicted format.

    :param filepath: the path to a file of unknown format
    :type filepath: pathlib.Path
    :return: fmt
    """
    fmt = None

    with open(filepath) as file_sniffer:
        line = file_sniffer.readline()

    # For now, sniffer only knows about FASTA and Genbank flatfile formats
    if line.startswith(">"):
        fmt = "fasta"
    elif line.startswith("LOCUS"):
        fmt = "genbank"

    return fmt
