from dna_features_viewer import BiopythonTranslator

from depht.classes.prophage import DEFAULT_PRODUCT

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FONT_FAMILY = "monospace"


class CircularSourceFeatureTranslator(BiopythonTranslator):
    """Translator for Circular Graphic."""

    def compute_feature_color(self, feature):
        """Compute the color of the feature."""
        if feature.location.strand == -1:   # reverse
            source_color = "red"
        else:   # forward
            source_color = "#449f66"

        color_map = {
            "source": source_color,  # light orange
            "gene": "red",
            "tRNA": "yellow",
            "CDS": "red",  # pink
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

    def compute_feature_fontdict(self, feature):
        """Compute a  font dict for this feature."""
        return {"family": DEFAULT_FONT_FAMILY}

    def compute_feature_box_linewidth(self, feature):
        """Compute a box_linewidth for this feature."""
        return 1.0


class LinearFeatureTranslator(BiopythonTranslator):
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
    att_types = ["attL", "attR"]
    ignored_features_types = ["source", "gene"]
    label_fields = ["product", "name", "note", "gene"]

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
            label = "|".join([str(x) for x in label])

        if label in self.hyp:
            label = None
        elif feature.type == "misc_recomb":
            for att_type in self.att_types:
                if att_type in label:
                    label = att_type
                    break

        return label

    def compute_feature_html(self, feature):
        """Compute the tooltip display text of the feature"""
        label = self.compute_feature_label(feature)
        properties = {"gb_type": feature.type}

        if feature.type == "CDS":
            gene = str(feature.qualifiers["gene"][0])
            translation = feature.qualifiers["translation"][0]
            if label is None:
                label = DEFAULT_PRODUCT

            properties["gene"] = gene
            properties["product"] = label
            properties["translation"] = translation

        elif feature.type == "misc_recomb":
            sequence = feature.qualifiers["note"][0]

            properties["name"] = label
            properties["sequence"] = sequence

        elif feature.type in ("tRNA", "tmRNA"):
            if label is None:
                label = ""

            properties["product"] = label
            properties["codon"] = feature.qualifiers.get("note", [""])[0]

        return properties

    def compute_feature_fontdict(self, feature):
        """Compute a  font dict for this feature."""
        return {"family": DEFAULT_FONT_FAMILY}

    def compute_feature_label_link_color(self, feature):
        """Compute the color of the  line linking the label to its feature"""
        return None

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
            CDS_color = "red"

        color_dict = {"source": "#34ff56",
                      "gene": "#f45562",
                      "CDS": CDS_color,
                      "misc_recomb": "#fbf3f6",
                      "tRNA": '#4554f4',
                      "tmRNA": '#56f443'}

        return color_dict[feature.type]
