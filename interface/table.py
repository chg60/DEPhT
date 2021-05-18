"""Create a table for the display."""

from tabulate import tabulate
from Bio import SeqIO

tabulate.PRESERVE_WHITESPACE = True


def get_features(file="output/GD43A_prophages.gb"):
    """
    Put features of a prophage feature file into a dictionary.

    :param file: file with prophage features
    :type file: str
    :return: characteristics of prophages
    :rtype: dict
    """
    info = {
        "Prophage Name\t": [],
        "Left Coordinate": [],
        "Right Coordinate": [],
        "Length": [],
        "Orientation": []}
    data = SeqIO.read(file, "gb")
    for feature in data.features:
        # add the name
        info.get("Prophage Name\t").append(
            feature.qualifiers.get("locus_tag")[0])
        # add the location Left Coordinates
        info.get("Left Coordinate").append(feature.location.nofuzzy_start)
        # add Location Right Coordinates
        info.get("Right Coordinate").append(feature.location.nofuzzy_end)
        # add length
        info.get("Length").append(
            feature.location.nofuzzy_end -
            feature.location.nofuzzy_start)
        # add orientation
        info.get("Orientation").append(feature.location.strand)
    return info


def make_table(data_dictionary):
    """
    Make a table.

    :param data_dictionary: file with prophage features
    :type data_dictionary: dict
    """
    table = tabulate(data_dictionary, headers='keys', tablefmt='html')
    # print(table)
    return table


def main():
    """Run the functions."""
    information = get_features()
    make_table(information)


main()
