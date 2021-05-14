"""Interactive diplay of genome maps for Prophicient."""

from bokeh.embed import file_html
from bokeh.resources import CDN
from pathlib import Path
import argparse

from linear_plot import CustomTranslator
import circular_display


def circular_graph_list(dir_path):
    """Get the filepaths for the svg files.

    :param dir_path: path to directory where genbank files are located
    :type dir_path: Path
    """
    dir_name = dir_path.stem

    circular_list = circular_display.get_path(dir_name)

    return circular_list


# translator = CustomTranslator()
def draw_figures(dir_path, filename, genome_list):
    """
    Plot the genome.

    :param dir_path: path to directory where svg file is located
    :type dir_path: Path
    :param filename: name of html file
    :type filename: str
    :param genome_name: name of the genome
    :type genome_name: str
    """
    circular_list = circular_graph_list(dir_path)

    translator = CustomTranslator()
    html_filename = filename + ".html"

    with open(html_filename, "w") as f:

        # write svg to html file
        for graph in circular_list:
            img_write = "<img src = \"" + str(graph)
            img_write += "\" style=\"height:750px\">"
            f.write(img_write)

        for genome in genome_list:
            record = translator.translate_record(str(genome))
            bokeh_plot = record.plot_with_bokeh(figure_width=100,
                                                figure_height="auto",
                                                tools="auto")
            html = file_html(bokeh_plot, CDN, filename)
            title = "<h2 style = \"font-family: monospace\">"
            title += str(genome.stem) + "</h2>"
            f.write(title)
            f.write(html)
            f.write("<br>")
            f.write("<br>")


def generate_prophage_list(dir_path):
    """
    Generate a list of genome names from the filepath.

    :parm dir_path: Path where genbank files are located
    :type dir_path: Path
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
    parser.add_argument("-f", "--filename", type=str,
                        default="results",
                        help=filename_help)
    args = parser.parse_args()

    if args.filename == "results":
        args.filename = args.dir.stem

    dir_path = args.dir/"prophages"

    genome_list = generate_prophage_list(dir_path)

    draw_figures(args.dir, args.filename, genome_list)


if __name__ == '__main__':
    main()
