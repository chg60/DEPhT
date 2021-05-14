"""Circular Display for Prophicient."""

from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
from pathlib import Path
import sys


class ExpressionUnitTranslator(BiopythonTranslator):
    """Translator for Circular Garphic."""

    def compute_feature_color(self, feature):
        """Compute the color of the feature."""
        if feature.location.strand == -1:
            color_map = {
                "source": "red",  # light orange
                "gene": "red",
                "tRNA": "yellow",
                "CDS": "red",  # pink
                "misc_recomb": "darkblue",
                "misc_feature": "#d1e9f1",  # light blue
                "backbone": "darkblue",
            }
            return color_map[feature.type]
        else:
            color_map = {
                "source": "#45f432",  # light orange
                "gene": "red",
                "tRNA": "yellow",
                "CDS": "#45f432",
                "misc_recomb": "darkblue",
                "misc_feature": "#d1e9f1",  # light blue
                "backbone": "darkblue",
            }
            return color_map[feature.type]

    def compute_feature_label(self, feature):
        """Compute the label of the feature."""
        if feature.type not in ("source"):
            return None
        else:
            return feature.qualifiers.get("locus_tag")[0]

    def compute_feature_box_linewidth(self, feature):
        """Compute a box_linewidth for this feature."""
        return 1.0

    # def compute_feature_legend_text(self, feature):
        # """Compute the legend text for the feature."""
        # return NAME


def get_path(folder_name):
    """Get the path of the folder.

    :param folder_name: Name of the folder
    :type folder_name: str
    """
    dir_name = Path.cwd()
    path_name = dir_name / folder_name

    figures = []

    for file in path_name.iterdir():
        # NAME = str(file.name)
        if file.suffix == ".gb":
            stem = str(file)
            translator = ExpressionUnitTranslator()
            graphic_record = translator.translate_record(
                stem, record_class=CircularGraphicRecord
            )
            ax, _ = graphic_record.plot(figure_width=10)
            figure_name = file.stem + "_circular_graphic.svg"
            figure_path = dir_name / figure_name
            ax.figure.savefig(figure_name, bbox_inches="tight")
            figures.append(figure_path)

    return figures


def main():
    """Run the functions."""
    get_path(sys.argv[1])


if __name__ == '__main__':
    main()
