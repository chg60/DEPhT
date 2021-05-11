from bokeh.embed import file_html
from bokeh.resources import CDN
from linear_plot import CustomTranslator

# translator = CustomTranslator()
def plot(genome_name):
    """
    Plot the genome.

    :param genome_name: name of the genome
    :type genome_name: str
    """
    translator = CustomTranslator()
    record = translator.translate_record(genome_name + ".gb")
    bokeh_plot = record.plot_with_bokeh(figure_width=50,
                                        figure_height=1,
                                        tools="auto")
    html = file_html(bokeh_plot, CDN, "my plot")

    html_filename = genome_name + ".html"

    with open(html_filename, "w") as f:
        f.write(html)


def main():
    """Run functions."""
    filename = input("enter name: ")
    plot(filename)


if __name__ == '__main__':
    main()
