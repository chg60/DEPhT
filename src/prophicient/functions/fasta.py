import pathlib


def parse_fasta(filepath):
    """
    Parses a fasta file into a list of headers and a corresponding list
    of sequences.

    :param filepath: fasta file to parse
    :type filepath: pathlib.Path or str
    :return: headers, sequences
    """
    headers, sequences = list(), list()

    fasta_reader = open(filepath, "r")

    cache = list()
    for line in fasta_reader:
        # If header line, flush the cache and store header
        if line.startswith(">"):
            sequences.append("".join(cache))
            cache = list()
            headers.append(line.lstrip(">").rstrip())
        # Otherwise, append to the cache
        else:
            cache.append(line.rstrip())

    # Flush the last sequence out of the cache, and pop empty sequence
    sequences.append("".join(cache))
    sequences.pop(0)

    # Close the file handle
    fasta_reader.close()

    return headers, sequences


def write_fasta(headers, sequences, filepath, width=80):
    """
    Writes the given headers and sequences to filepath.

    :param headers: the sequence labels
    :type headers: list of str
    :param sequences: the sequences to write
    :type sequences: list of str
    :param filepath: the fasta file to write
    :type filepath: pathlib.Path or str
    :param width: maximum number of characters per sequence line
    :type width: int
    """
    if not isinstance(headers, list) or not isinstance(sequences, list):
        raise TypeError(f"headers and sequences should be lists of strings")

    fasta_writer = open(filepath, "w")

    for header, sequence in zip(headers, sequences):
        fasta_writer.write(f">{header}\n")
        # Use string slicing to split the sequence to satisfy width param
        for i in range(0, len(sequence), width):
            fasta_writer.write(f"{sequence[i:i + width]}\n")

    # Close the file handle
    fasta_writer.close()
