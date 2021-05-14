"""BiopythonTrasnlator modified for use for Prophicient display."""

from dna_features_viewer import BiopythonTranslator
from dna_features_viewer import GraphicFeature


class CustomTranslator(BiopythonTranslator):
    """Custom Translator adapted from BiopythonTranslator."""

    # Label fields indicates the order in which annotations fields are
    # considered to determine the feature's label

    # MODIFICATIONS:
    # label based on feature type!
    # display CDS unless HP, then --> gene
    # add display for tRNA, tmRNA - display note field
    # misc_recomb -- feature type - show sequence

    # bird's eye view of entire genome - circular
    # display - source?

    # feature.location - +/- for forward or reverse
    # feature.location.strand

    hyp = ["Hypothetical Protein", "hypothetical protein"]
    ignored_features_types = ["source", "gene"]
    label_fields = ["product", "note", "gene"]

    def compute_feature_box_linewidth(self, feature):
        """Compute a box_linewidth for this feature."""
        return 1.0

    def compute_feature_linewidth(self, feature):
        """Compute the edge width of the feature's arrow/rectangle."""
        return 10

    def compute_feature_label(self, feature):
        """Compute the label of the feature."""
        label = feature.type
        for key in self.label_fields:
            if key in feature.qualifiers and len(feature.qualifiers[key]):
                label = feature.qualifiers[key]
                break

        if isinstance(label, list):
            label = "|".join(label)

        if label in self.hyp:
            label = None

        return label

    def compute_feature_legend_text(self, feature):
        """Compute the feature legend text."""
        return feature.type

    def compute_feature_box_color(self, feature):
        """Compute a box_color for this feature."""
        return "black"

    def compute_feature_color(self, feature):
        """Compute color for the given feature - differences based on type."""
        if feature.location.strand == 1:
            # green for forward
            CDS_color = "#449f66"
        else:
            # red for reverse
            CDS_color = "#ff4466"

        color_dict = {"source": "#34ff56",
                      "gene": "#f45562",
                      "CDS": CDS_color,
                      "misc_recomb": "#fbf3f6",
                      "tRNA": '#4554f4',
                      "tmRNA": '#56f443'}

        return color_dict[feature.type]
