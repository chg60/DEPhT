import pathlib


def parse_fasta(f):
    """
    Parses a fasta file into a list of geneids and a corresponding
    list of translations.

    :param f: fasta file to parse
    :type f: pathlib.Path
    :return: geneids, translations
    """
    geneids = list()
    translations = list()
    with f.open("r") as fh:
        line = fh.readline()
        while line:
            if line.startswith(">"):
                geneids.append(line.lstrip(">").rstrip())
                translations.append("")
            else:
                translations[-1] += line.rstrip()
            line = fh.readline()

    return geneids, translations


def write_fasta(f, geneids, translations):
    """
    Writes the given geneids and translations to file f.

    Round trip with `write_fasta` and `parse_fasta` will produce the same result.

    :param f: the fasta file to write
    :type f: pathlib.Path
    :param geneids: the gene labels
    :type geneids: list
    :param translations: the gene sequences
    :type translations: list
    :return:

    >>> gs = ["Muddy_gp28", "Muddy_gp29"]
    >>> ts = ["MTKFWEFVKSDKFRLYFYSVCVAVMGLLVYYGIVEAEAVPYWLTLLGAIGMVGNATAAANLGSQIRQKGGEG",
    ... "MTWWSADFWNNIGPVGLSVLACVFFVVALVRGWLVIGRYHRETVERLDARAQKDAETIDVLSRAVTEKVAEDQATTRILSAIRDLWTSSKEEAT"]
    >>> fasta_file = pathlib.Path("/tmp/test.fasta")
    >>> write_fasta(fasta_file, gs, ts)
    >>> read_gs, read_ts = parse_fasta(fasta_file)
    >>> assert gs[0] == read_gs[0]
    >>> assert gs[1] == read_gs[1]
    >>> assert ts[0] == read_ts[0]
    >>> assert ts[1] == read_ts[1]
    """
    with f.open("w") as fh:
        for geneid, translation in zip(geneids, translations):
            fh.write(f">{geneid}\n{translation}\n")
