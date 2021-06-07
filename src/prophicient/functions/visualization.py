from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from bokeh.embed import file_html
from bokeh.resources import CDN
from dna_features_viewer import CircularGraphicRecord
from reportlab.lib import colors
from tabulate import tabulate

from prophicient.classes.file_translator import (
                                        CircularSourceFeatureTranslator,
                                        LinearFeatureTranslator)

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_NAME = "cartoon"

# PRESETS
# -----------------------------------------------------------------------------
tabulate.PRESERVE_WHITESPACE = True


# BOKEH/DNAFEATUREVIEWER FUNCTIONS
# -----------------------------------------------------------------------------
def draw_complete_diagram(outdir, contigs, prophage_records, tmp_dir,
                          name=DEFAULT_NAME):
    """
    Plots a circular genome diagram of the host bacterial record(s), linear
    genome diagrams for the predicted prophage record(s), and a table
    containing metadata assocated with their relationship.

    :param outdir: File directory the graphics will be written to
    :type outdir: pathib.Path
    :type contigs: Auto-annotated bacterial sequence contigs
    :type contigs: list
    :param prophage_paths: Filepaths of the associated prophages
    :type prophage_paths: list
    :param tmp_dir: Path to place modified file resources
    :type tmp_dir: pathib.Path
    """
    # Scrub each bacterial record contig of non-source features
    # and write to file
    host_record_paths = scrub_host_records(contigs, tmp_dir)

    resources_dir = outdir.joinpath("resources")
    resources_dir.mkdir(exist_ok=True)

    # Scrape metadata from each genome contig record
    host_metadata = scrape_and_tabulate_host_metadata(host_record_paths)

    # Create a circular image resource for each genome contig record
    host_graphic_resources = create_host_graphic_resources(contigs,
                                                           resources_dir)

    # Embed prophage plots
    embedded_prophage_plots = embed_prophage_genomes(prophage_records)

    filepath = outdir.joinpath(name).with_suffix(".html")

    with filepath.open(mode="w") as filehandle:
        for i in range(len(contigs)):
            graph = host_graphic_resources[i]
            table = host_metadata[i]

            img_html = ("".join(["<img src = \"",
                                 str(graph.relative_to(outdir)),
                                 "\"style=\"height:750px\">"]))
            filehandle.write(img_html)

            table_html = ("".join(["<br>", table, "<br>"]))
            filehandle.write(table_html)

        for i in range(len(prophage_records)):
            title = str(prophage_records[i].id)

            title_html = "".join(["<h2_style = \"font-family: monospace\">",
                                  title, "</h2>"])
            filehandle.write(title_html)

            plot = embedded_prophage_plots[i]
            filehandle.write(plot)

            filehandle.write("<br><br>")


def scrub_host_records(contigs, tmp_dir):
    host_record_paths = []
    for contig in contigs:
        source_features = []

        # Scrape all source features
        for feature in contig.features:
            if feature.type != "source":
                continue

            # Ignore source features describing the whole contig
            if len(feature) >= len(contig.seq):
                continue

            source_features.append(feature)

        contig.features = source_features

        record_path = tmp_dir.joinpath(str(contig.id)).with_suffix(".gbk")
        host_record_paths.append(record_path)

        SeqIO.write(contig, record_path, "genbank")

    return host_record_paths


def scrape_and_tabulate_host_metadata(host_record_paths):
    """
    Put features of a prophage feature files into dictionaries.

    :param host_record_paths: list of files with prophage features
    :type file: list
    :return: characteristics of prophages
    :rtype: list
    """
    host_metadata = []
    for host_record_path in host_record_paths:
        info = {
            "Prophage Name\t": [],
            "Left Coordinate": [],
            "Right Coordinate": [],
            "Length": [],
            "Orientation": []}
        data = SeqIO.read(host_record_path, "gb")
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

        table = tabulate(info, headers="keys", tablefmt="html")
        host_metadata.append(table)

    return host_metadata


def create_host_graphic_resources(contigs, resources_dir):
    host_graphic_resources = []
    for contig in contigs:
        host_translator = CircularSourceFeatureTranslator()

        graphic_record = host_translator.translate_record(
                                            contig,
                                            record_class=CircularGraphicRecord)

        ax, _ = graphic_record.plot(figure_width=10)

        resource_filepath = resources_dir.joinpath(contig.id).with_suffix(
                                                                        ".svg")

        ax.figure.savefig(resource_filepath, bbox_inches="tight")
        host_graphic_resources.append(resource_filepath)

    return host_graphic_resources


def embed_prophage_genomes(prophage_records):
    translator = LinearFeatureTranslator()

    embedded_prophage_plots = []
    for record in prophage_records:
        prophage_translation = translator.translate_record(record)
        bokeh_plot = prophage_translation.plot_with_bokeh(figure_width=100,
                                                          figure_height="auto",
                                                          tools="auto")

        embedded_plot = file_html(bokeh_plot, CDN, record.id)
        embedded_prophage_plots.append(embedded_plot)

    return embedded_prophage_plots


# REPORTLAB FUNCTIONS
# -----------------------------------------------------------------------------
def prophage_diagram(record, filepath):
    """
    Plots a linear genome diagram of the given record, and saves it
    in the indicated file.

    :param record: the record to plot
    :type record: Bio.SeqRecord.SeqRecord
    :param filepath: the file to store the genome diagram in
    :type filepath: pathlib.Path
    """
    diagram = GenomeDiagram.Diagram(f"{record.id}")
    track = diagram.new_track(track_level=1, name=f"{filepath.stem}",
                              greytrack=True, greytrack_labels=1,
                              scale_smalltick_interval=1000,
                              scale_largetick_interval=5000)
    track_features = track.new_set()

    for feature in record.features:
        add_track_feature(track_features, feature)

    diagram.draw(format="linear", orientation="landscape",
                 pagesize="A4", fragments=8, start=0, end=len(record)-1)
    diagram.write(str(filepath), "PDF")


def add_track_feature(feature_set, feature):
    if feature.type not in ("CDS", "tRNA", "tmRNA"):
        return

    if feature.type == "CDS":
        sigil = "ARROW"
        name = feature.qualifiers["product"][0]
        if name == "hypothetical protein":
            name = ""
    elif feature.type == "tRNA":
        sigil, color = "BOX", colors.blue
        name = feature.qualifiers["note"][0]
    else:
        sigil, color = "BOX", colors.gold
        name = feature.type

    if feature.location.strand == 1:
        color, angle = colors.green, 45
    else:
        color, angle = colors.red, 225

    feature_set.add_feature(feature, color=color, label=True, name=name,
                            labelsize=3, label_angle=angle, sigil=sigil,
                            label_position="middle")
