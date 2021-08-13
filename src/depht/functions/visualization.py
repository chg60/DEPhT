import pandas
import pretty_html_table
import bokeh
from bokeh import plotting
from bokeh.embed import file_html
from bokeh.models import HoverTool, Range1d, Title
from bokeh.resources import CDN
from dna_features_viewer import CircularGraphicRecord
from matplotlib import pyplot

from depht.classes.file_translator import (
    CircularSourceFeatureTranslator, LinearFeatureTranslator,
    DEFAULT_FONT_FAMILY)

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_NAME = "cartoon"


# BOKEH/DNAFEATUREVIEWER FUNCTIONS
# -----------------------------------------------------------------------------
def draw_complete_diagram(outdir, contigs, prophages, tmp_dir,
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
    scrub_host_records(contigs)

    # Scrape metadata from each genome contig record
    host_metadata = scrape_and_tabulate_host_metadata(contigs)

    # Create a circular image resource for each genome contig record
    host_graphic_resources = create_host_graphic_resources(contigs,
                                                           tmp_dir)

    contig_prophage_map = {}
    for prophage in prophages:
        contig_prophages = contig_prophage_map.get(prophage.parent_record.id,
                                                   list())
        contig_prophages.append(prophage)
        contig_prophage_map[prophage.parent_record.id] = contig_prophages

    filepath = outdir.joinpath(f"{name}.html")

    filehandle = open(filepath, "w")
    for i in range(len(contigs)):
        contig = contigs[i]
        if not contig.features:
            continue

        graph = host_graphic_resources[i]
        table = host_metadata[i]

        filehandle.write(graph)

        table_html = ("".join([table, "<br>"]))
        filehandle.write(table_html)

        contig_prophages = contig_prophage_map.get(contig.id, list())

        # Embed prophage plots
        embedded_prophage_plots = embed_prophage_genomes(contig_prophages)

        for prophage_plot in embedded_prophage_plots:
            filehandle.write(prophage_plot)
            filehandle.write("<br><br>")

    filehandle.close()

    bokeh.io.reset_output()


def scrub_host_records(contigs):
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


def scrape_and_tabulate_host_metadata(contigs):
    """
    Put features of a prophage feature files into dictionaries.

    :param host_record_paths: list of files with prophage features
    :type file: list
    :return: characteristics of prophages
    :rtype: list
    """
    host_metadata = []
    for contig in contigs:
        keys = ["Prophage Name", "Left Coordinate", "Right Coordinate",
                "Length"]
        fwd_info = {key: list() for key in keys}
        rvs_info = {key: list() for key in keys}

        for feature in contig.features:
            if feature.location.strand == -1:
                info = rvs_info
            else:
                info = fwd_info
            # add the name
            info.get("Prophage Name").append(
                feature.qualifiers.get("locus_tag")[0])
            # add the location Left Coordinates
            info.get("Left Coordinate").append(feature.location.nofuzzy_start)
            # add Location Right Coordinates
            info.get("Right Coordinate").append(feature.location.nofuzzy_end)
            # add length
            info.get("Length").append(
                feature.location.nofuzzy_end -
                feature.location.nofuzzy_start)

        fwd_table = pandas.DataFrame(fwd_info)
        rvs_table = pandas.DataFrame(rvs_info)

        if len(fwd_table) > 0:
            fwd_table_html = pretty_html_table.build_table(
                                            fwd_table, "green_light",
                                            font_size="14px",
                                            text_align="center",
                                            font_family=DEFAULT_FONT_FAMILY)
        else:
            fwd_table_html = ""

        if len(rvs_table) > 0:
            rvs_table_html = pretty_html_table.build_table(
                                            rvs_table, "red_light",
                                            font_size="14px",
                                            text_align="center",
                                            font_family=DEFAULT_FONT_FAMILY)
        else:
            rvs_table_html = ""

        table_html = "\n".join([fwd_table_html, rvs_table_html])
        host_metadata.append(table_html)

    return host_metadata


def create_host_graphic_resources(contigs, resources_dir):
    host_graphic_resources = []
    for contig in contigs:
        host_translator = CircularSourceFeatureTranslator()

        graphic_record = host_translator.translate_record(
                                            contig,
                                            record_class=CircularGraphicRecord)

        ax, _ = graphic_record.plot(figure_width=10, with_ruler=True)

        # Write center label
        ax.text(0, -1, contig.id,
                horizontalalignment="center", verticalalignment="center",
                fontfamily="monospace", fontsize=30, fontweight="demibold")

        # Write center bp
        ax.text(0, -1.1,
                " ".join([str(len(contig)), "bp"]),
                horizontalalignment="center", verticalalignment="center",
                fontfamily="monospace", fontsize=12, fontweight="normal")

        resource_filepath = resources_dir.joinpath(contig.id).with_suffix(
                                                                        ".svg")
        ax.figure.savefig(resource_filepath, bbox_inches="tight")
        pyplot.close(ax.figure)

        filehandle = open(resource_filepath, "r")

        lines = filehandle.readlines()
        lines = "".join(lines)
        host_graphic_resources.append(lines)

        filehandle.close()

    return host_graphic_resources


def embed_prophage_genomes(prophages):
    translator = LinearFeatureTranslator()

    embedded_prophage_plots = []
    for prophage in prophages:
        graphic_record = translator.translate_record(prophage.record)
        bokeh_plot = plot_with_bokeh(graphic_record, figure_width=100,
                                     figure_height="auto", tools="auto",
                                     label_link_width=0,
                                     label_text_font=DEFAULT_FONT_FAMILY)

        # Format primary title
        bokeh_plot.title.text = prophage.id
        bokeh_plot.title.align = "left"
        bokeh_plot.title.text_font = DEFAULT_FONT_FAMILY
        bokeh_plot.title.text_font_size = "18px"

        # Add secondary location title
        location_prefix = ""
        if prophage.strand == -1:
            location_prefix = "complement"

        location = "".join([location_prefix, "(", str(prophage.start),
                            "..", str(prophage.end), ")"])
        bokeh_plot.add_layout(Title(text=location, align="left",
                                    text_font_size="12px",
                                    text_font=DEFAULT_FONT_FAMILY), "above")

        embedded_plot = file_html(bokeh_plot, CDN, prophage.id)
        embedded_prophage_plots.append(embedded_plot)

    return embedded_prophage_plots


# BOKEH HELPER OVERRIDE
# -----------------------------------------------------------------------------

# If you're wondering why there's a redundant function here,
# ask why dna_features_viewer.GraphicRecord.plot_with_bokeh hard codes so much
def plot_with_bokeh(graphic_record, figure_width=100, figure_height="auto",
                    title=None, tools="auto", feature_box_color="#000000",
                    label_link_color="#000000", label_text_font="arial",
                    label_text_font_size="12px",
                    label_text_font_style="normal", label_link_width=0.5):
    """Plots the graphic record using Bokeh"""
    if tools == "auto":
        tools = "xpan,xwheel_zoom,reset,tap"

    # Scrape matplotlib plot info

    ax, (features_levels, plot_data) = graphic_record.plot(
                                                figure_width=figure_width,
                                                annotate_inline=False)
    width, height = [int(100 * e) for e in ax.figure.get_size_inches()]
    pyplot.close(ax.figure)

    if figure_height == "auto":
        height = int(0.5 * height)
    else:
        height = 100 * figure_height

    # Ensures that height of the plot is 185 or greater, the minimal height
    # to see all icons
    height = max(height, 185)

    max_y = max(
                [data["annotation_y"] for f, data in plot_data.items()]
                + list(features_levels.values()))

    # Build the Bokeh plot
    plot = plotting.figure(plot_width=width, plot_height=height, tools=tools,
                           x_range=Range1d(0, graphic_record.sequence_length),
                           y_range=Range1d(-1, max_y + 1))

    # Define different hovertools
    cds_tooltips = [("gene", "@gene"), ("product", "@product"),
                    ("location", "@location")]
    cds_hovertool = HoverTool(tooltips=cds_tooltips, mode="mouse",
                              names=["cds"])
    plot.add_tools(cds_hovertool)

    att_tooltips = [("name", "@name"), ("sequence", "@sequence")]
    att_hovertool = HoverTool(tooltips=att_tooltips, mode="mouse",
                              names=["att"])
    plot.add_tools(att_hovertool)

    trna_tooltips = [("product", "@product"), ("codon", "@codon"),
                     ("location", "@location")]
    trna_hovertool = HoverTool(tooltips=trna_tooltips, mode="mouse",
                               names=["tRNA", "tmRNA"])
    plot.add_tools(trna_hovertool)

    # Set up patches source
    cds_feature_patches = list()
    att_feature_patches = list()
    trna_feature_patches = list()
    for feature, level in features_levels.items():
        bokeh_feature_patch = graphic_record.bokeh_feature_patch(
                                    feature.start, feature.end, feature.strand,
                                    figure_width=figure_width, level=level,
                                    color=feature.color, label=feature.label)

        location_prefix = ""
        if feature.strand == -1:
            location_prefix = "complement"

        location = "".join([location_prefix, "(", str(feature.start),
                            "..", str(feature.end), ")"])

        html_data = feature.html
        if html_data["gb_type"] == "CDS":
            bokeh_feature_patch["gene"] = html_data["gene"]
            bokeh_feature_patch["product"] = html_data["product"]
            bokeh_feature_patch["location"] = location
            cds_feature_patches.append(bokeh_feature_patch)
        elif feature.html["gb_type"] == "misc_recomb":
            bokeh_feature_patch["name"] = html_data["name"]
            bokeh_feature_patch["sequence"] = html_data["sequence"]
            att_feature_patches.append(bokeh_feature_patch)
        elif feature.html["gb_type"] in ("tRNA", "tmRNA"):
            bokeh_feature_patch["product"] = html_data["product"]
            bokeh_feature_patch["codon"] = html_data["codon"]
            bokeh_feature_patch["location"] = location
            trna_feature_patches.append(bokeh_feature_patch)

    if cds_feature_patches:
        cds_patches_source = plotting.ColumnDataSource(
                                            pandas.DataFrame.from_records(
                                                    cds_feature_patches))
        plot.patches(xs="xs", ys="ys", color="color", name="cds",
                     line_color=feature_box_color, source=cds_patches_source)

    if att_feature_patches:
        att_patches_source = plotting.ColumnDataSource(
                                            pandas.DataFrame.from_records(
                                                    att_feature_patches))
        plot.patches(xs="xs", ys="ys", color="color", name="att",
                     line_color=feature_box_color, source=att_patches_source)

    if trna_feature_patches:
        trna_patches_source = plotting.ColumnDataSource(
                                            pandas.DataFrame.from_records(
                                                    trna_feature_patches))
        plot.patches(xs="xs", ys="ys", color="color", name="tRNA",
                     line_color=feature_box_color, source=trna_patches_source)

    # Set up floating label and connector source
    bokeh_floating_labels = list()
    bokeh_label_segments = list()
    for feature, pdata in plot_data.items():
        bokeh_floating_label = dict(
                                x=feature.x_center, y=pdata["annotation_y"],
                                text=feature.label, color=feature.color)
        bokeh_floating_labels.append(bokeh_floating_label)

        bokeh_segment = dict(x0=feature.x_center, x1=feature.x_center,
                             y0=pdata["annotation_y"], y1=pdata["feature_y"])
        bokeh_label_segments.append(bokeh_segment)

    if plot_data != {}:
        labels_source = plotting.ColumnDataSource(
                                        pandas.DataFrame.from_records(
                                                    bokeh_floating_labels))
        plot.text(x="x", y="y", text="text", text_align="center",
                  text_font_size=label_text_font_size,
                  text_font=label_text_font,
                  text_font_style=label_text_font_style,
                  source=labels_source)

        segments_source = plotting.ColumnDataSource(
                                        pandas.DataFrame.from_records(
                                                    bokeh_label_segments))
        plot.segment(x0="x0", x1="x1", y0="y0", y1="y1",
                     line_width=label_link_width, color=label_link_color,
                     source=segments_source)

    plot.yaxis.visible = False
    plot.outline_line_color = None
    plot.grid.grid_line_color = None
    plot.toolbar.logo = None

    return plot
