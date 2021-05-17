"""Interactive display of genome maps for Prophicient."""

from bokeh.embed import file_html
from bokeh.resources import CDN
from pathlib import Path
from linear_plot import CustomTranslator
import circular_display
import path_to_directories
import table


def circular_graph_list(dir_path):
    """Get the filepaths for the svg files.

    :param dir_path: path to directory where genbank files are located
    :type dir_path: Path
    """
    dir_name = dir_path.stem

    circular_list = circular_display.get_path(dir_name)

    return circular_list


def draw_circular(dir_list):
    """Draw the circular graphs.

    :param dir_list: List of directories
    :type dir_list: list
    """
    circular_list = circular_display.get_path(dir_list)

    return circular_list


def draw_linear(translator, prophage):
    """Draw the linear graphics.

    :param translator: User-definde BiopythonTranslator
    :type translator: CustomTranslator
    :param prophage: path to prophage
    :type prophage: Path
    :returns bokeh_plot: plotted bokeh plot
    :rtype bokeh_plot: bokeh.plotting.figure.Figure
    """
    record = translator.translate_record(str(prophage))
    bokeh_plot = record.plot_with_bokeh(figure_width=100,
                                        figure_height="auto",
                                        tools="auto")
    return bokeh_plot


# translator = CustomTranslator()
def draw_figures(filename, dir_list):
    """
    Plot the genome.

    :param dir_path: path to directory where svg file is located
    :type dir_path: Path
    :param filename: name of html file
    :type filename: str
    :param genome_name: name of the genome
    :type genome_name: str
    """
    circular_list = draw_circular(dir_list[1])

    prophage_list = dir_list[0]

    table_info = table.get_features()
    proph_table = table.make_table(table_info)

    translator = CustomTranslator()
    html_filename = filename + ".html"

    with open(html_filename, "w") as f:
        # write svg to html file
        for graph in circular_list:
            img_write = "<img src = \"" + str(graph)
            img_write += "\" style=\"height:750px\">"
            f.write(img_write)

        f.write("<br>")
        f.write(proph_table)
        f.write("<br>")

        for proph in prophage_list:
            proph = Path(proph)
            bokeh_plot = draw_linear(translator, proph)
            html = file_html(bokeh_plot, CDN, filename)
            title = "<h2 style = \"font-family: monospace\">"
            title += str(proph.stem) + "</h2>"
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
    dir_list = path_to_directories.extract_files("prophages", "output")
    filename = "results"

    draw_figures(filename, dir_list)


if __name__ == '__main__':
    main()
