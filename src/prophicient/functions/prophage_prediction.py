import pickle

import pandas as pd
from Bio.SeqUtils import GC123
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from prophicient.__main__ import PACKAGE_DIR

MODEL_PATH = PACKAGE_DIR.joinpath("data/prophage_model.pickle")
RBS_SCORES = {"None": {"None": 0},
              "GGA/GAG/AGG": {"3-4bp": 1, "11-12bp": 6, "5-10bp": 13},
              "3Base/5BMM": {"13-15bp": 2},
              "4Base/6BMM": {"13-15bp": 3},
              "AGxAG": {"11-12bp": 4, "3-4bp": 5, "5-10bp": 9},
              "GGxGG": {"11-12bp": 4, "3-4bp": 8, "5-10bp": 14},
              "AGGAG(G)/GGAGG": {"13-15bp": 10},
              "AGGA/GGAG/GAGG": {"3-4bp": 11, "11-12bp": 12},
              "AGGA": {"5-10bp": 15},
              "GGAG/GAGG": {"5-10bp": 16},
              "AGxAGG/AGGxGG": {"11-12bp": 17, "3-4bp": 18, "5-10bp": 19},
              "AGGAG/GGAGG": {"11-12bp": 20},
              "AGGAG": {"3-4bp": 21, "5-10bp": 22},
              "GGAGG": {"3-4bp": 23, "5-10bp": 24},
              "AGGAGG": {"11-12bp": 25, "3-4bp": 26, "5-10bp": 27}}
FEATURE_KEYS = ("left", "right", "strand", "length", "upstream", "downstream",
                "gc1", "gc2", "gc3", "gravy", "isoelectric", "aromaticity")
# "rbs_score" not at all predictive, in Mycobacteria
# "C", "H", "Q" are most independent amino acid frequencies
# "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
# "S", "T", "V", "W", and "Y" are slightly correlated with other features

BACTERIA = 0        # Gene prediction state - bacterial
PROPHAGE = 1        # Gene prediction state - prophage


def contig_to_dataframe(contig):
    """
    Walks the features of the given contig and calculates the feature
    data to be used for model training and/or model prediction.

    :param contig: a contig to analyze features from
    :type contig: Bio.SeqRecord.SeqRecord
    :return: dataframe
    """
    gene_dict = dict()
    for key in FEATURE_KEYS:
        gene_dict[key] = list()

    cds_features = [feat for feat in contig.features if feat.type == "CDS"]
    for j, feature in enumerate(cds_features):
        left, right = feature.location.start, feature.location.end
        strand = feature.location.strand

        # Can't calculate first gene's upstream distance yet...
        if j == 0:
            upstream = 0
        else:
            upstream = left - gene_dict["right"][-1] + 1

        _, gc1, gc2, gc3 = GC123(str(feature.extract(contig).seq))

        # note = feature.qualifiers["note"][0].split("; ")
        # rbs_motif = note[0].split(": ")[-1]
        # rbs_spacer = note[1].split(": ")[-1]
        # rbs_score = RBS_SCORES[rbs_motif][rbs_spacer]

        p_analysis = ProteinAnalysis(feature.qualifiers["translation"][0])

        gene_dict["left"].append(left)
        gene_dict["right"].append(right)
        gene_dict["strand"].append(strand)
        gene_dict["length"].append(len(feature))
        gene_dict["upstream"].append(upstream)
        # gene_dict["rbs_score"].append(rbs_score)
        gene_dict["gc1"].append(gc1)
        gene_dict["gc2"].append(gc2)
        gene_dict["gc3"].append(gc3)
        gene_dict["gravy"].append(p_analysis.gravy())
        gene_dict["isoelectric"].append(p_analysis.isoelectric_point())
        gene_dict["aromaticity"].append(p_analysis.aromaticity())
        # for aa, freq in p_analysis.get_amino_acids_percent().items():
        #     gene_dict[aa].append(freq)

    # Fix gene_dict["upstream"][0]
    upstream = gene_dict["left"][0] - gene_dict["right"][-1] + 1 + len(contig)
    gene_dict["upstream"][0] = upstream

    # Create gene_dict["downstream"]
    downstream = [x for x in gene_dict["upstream"]]
    downstream.append(downstream.pop(0))
    gene_dict["downstream"] = downstream

    # Update gene_dict["strand"] to be something meaningful
    window_size = 10  # Local window size
    num_genes = 2 * window_size + 1
    strands = list()
    for i in range(len(gene_dict["strand"])):
        cursor_strand = None
        local_strands = []
        min_i, max_i = i - window_size, i + window_size + 1
        for j in range(min_i, max_i):
            if j < 0:
                j += len(gene_dict["strand"])
            elif j > len(gene_dict["strand"]) - 1:
                j -= len(gene_dict["strand"])
            local_strands.append(gene_dict["strand"][j])
            if i == j:
                cursor_strand = local_strands[-1]
        strand = float(local_strands.count(cursor_strand)) / num_genes
        strands.append(strand)
    gene_dict["strand"] = strands

    # Absolute left and right coordinates don't matter anymore
    gene_dict.pop("left"), gene_dict.pop("right")

    return pd.DataFrame(gene_dict)


def predict_prophage_genes(dataframe, model_path=MODEL_PATH):
    """
    Uses an sklearn classifier from `model_path` to predict the nature
    of genes contained in `dataframe`.

    A prediction will be made for each gene - 0 is a bacterial gene,
    1 is a prophage gene.

    :param dataframe: genes to make predictions for
    :type dataframe: pd.DataFrame
    :param model_path: path to a binary file with sklearn model inside
    :type model_path: pathlib.Path
    :return: predictions
    """
    with model_path.open("rb") as model_reader:
        model = pickle.load(model_reader)

    initial_predictions = [x[1] > 0.5 for x in model.predict_proba(dataframe)]

    # Smooth the predictions to enforce local continuity
    window_size = 25  # Local window size
    num_genes = 2 * window_size + 1
    predictions = list()
    for i in range(len(initial_predictions)):
        local_predictions = list()
        min_i, max_i = i - window_size, i + window_size + 1
        for j in range(min_i, max_i):
            if j < 0:
                j += len(initial_predictions)
            elif j > len(initial_predictions) - 1:
                j -= len(initial_predictions)
            local_predictions.append(initial_predictions[j])
        if sum(local_predictions) * 2 > num_genes:
            prediction = 1
        else:
            prediction = 0
        predictions.append(prediction)

    return predictions


def predict_prophage_coords(contig, predictions, extend_by=10000):
    """
    Uses gene status predictions (0/1 = bacterial/prophage) to estimate
    coordinates for any prophages that may be on the indicated contig.

    :param contig: the contig to scan for possible prophages
    :type contig: Bio.SeqRecord.SeqRecord
    :param predictions: the binary predictions for this contig's genes
    :type predictions: list
    :param extend_by: number of basepairs to extend the prophages by
    :type extend_by: int
    :return: prophage_coords
    """
    prophage_coords = list()
    left, right = None, None
    n_inflections = 0
    previous_state = predictions[-1]

    # predictions are only made on the cds features
    cds_features = [feat for feat in contig.features if feat.type == "CDS"]

    for i, current_state in enumerate(predictions):
        if previous_state == BACTERIA and current_state == PROPHAGE:
            left = cds_features[i].location.start - extend_by
            n_inflections += 1
        elif previous_state == PROPHAGE and current_state == BACTERIA:
            right = cds_features[i - 1].location.end + extend_by
            prophage_coords.append((left, right))
            left, right = None, None
            n_inflections += 1
        previous_state = current_state
    # Record the last one found
    if left or right:
        prophage_coords.append((left, right))

    # Deal with edge cases (4 major possibilities):
    # Case 0: whole contig is a prophage - last gene processed was a predicted
    # prophage gene, and we had no PROPHAGE <-> BACTERIAL inflections
    if previous_state == PROPHAGE and n_inflections == 0:
        left, right = 0, len(contig)
        prophage_coords.append((left, right))

    # Case 1: there's a prophage that falls off the left edge of a contig -
    # len(prophage_coords) >= 1, and no left coord for first coordinate tuple
    if len(prophage_coords) >= 1 and not prophage_coords[0][0]:
        left, right = 0, prophage_coords[0][-1]
        prophage_coords[0] = (left, right)

    # Case 2: there's a prophage that falls off the right edge of a contig -
    # len(prophage_coords) >= 1, and no right coord for last coordinate tuple
    if len(prophage_coords) >= 1 and not prophage_coords[-1][-1]:
        left, right = prophage_coords[-1][0], len(contig) - 1
        prophage_coords[-1] = (left, right)

    return prophage_coords
