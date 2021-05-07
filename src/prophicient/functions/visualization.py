from reportlab.lib import colors
from Bio.Graphics import GenomeDiagram


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
