"""Interactive diplay of genome maps for Prophicient."""

from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.layouts import column
from linear_plot import CustomTranslator
from pathlib import Path
import argparse


# translator = CustomTranslator()
def plot(filename, genome_list):
    """
    Plot the genome.

    :param genome_name: name of the genome
    :type genome_name: str
    """
    translator = CustomTranslator()
    html_filename = filename + ".html"
    with open(html_filename, "w") as f:
        for genome in genome_list:
            record = translator.translate_record(str(genome))
            bokeh_plot = record.plot_with_bokeh(figure_width=50,
                                                figure_height=1,
                                                tools="auto")
            html = file_html(bokeh_plot, CDN, filename)
            title = "<h2 style = \"font-family: monospace\">"
            title += str(genome.stem) + "</h2>"
            f.write(title)
            f.write(html)
            f.write("<br>")
            f.write("<br>")


def generate_genome_list(dir_path):
    """
    Generate a list of genome names from the filepath.

    :parm filepath: Path where genbank files are located
    :type filepath: Path
    :returns: list of genome filepaths
    :rtype: list
    """
    genomes = []

    # make a list of all the genbank files
    for file in sorted(dir_path.iterdir()):
        if file.suffix == ".gb":
            genomes.append(file)

    return genomes


def main():
    """Run functions."""
    parser = argparse.ArgumentParser()
    dir_help = "Directory with genbank files"
    filename_help = "Name of file to store results to"

    parser.add_argument("dir", type=Path, help=dir_help)
    parser.add_argument("-f", "--filename", type=str, default="results",
                        help=filename_help)
    args = parser.parse_args()

    genome_list = generate_genome_list(args.dir)

    plot(args.filename, genome_list)


if __name__ == '__main__':
    main()
