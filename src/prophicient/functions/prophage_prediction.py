import pickle
import pandas as pd

from prophicient import PACKAGE_DIR
from prophicient.functions.sliding_window import *
from prophicient.functions.statistics import average

MODEL_PATH = PACKAGE_DIR.joinpath("data/prophage_model.pickle")

WINDOW = 51         # Number of CDS features to consider in a window

BACTERIA = 0        # Gene prediction state - bacterial
PROPHAGE = 1        # Gene prediction state - prophage


def average_gene_size(starts, stops, length, window=WINDOW):
    """
    Calculates and returns the average gene size in leading/centered/
    lagging windows of genes.

    :param starts: gene left coordinates
    :type starts: list of int
    :param stops: gene right coordinates
    :type stops: list of int
    :param length: contig length in nucleotides
    :type length: int
    :param window: how many genes to average over
    :type window: int
    :return: leading_sizes, center_sizes, lagging_sizes
    """
    if not len(starts) == len(stops):
        raise ValueError(f"len(starts) ({len(starts)}) != len(stops) ("
                         f"{len(stops)})")

    num_genes = len(starts)

    leading_sizes = list()
    lagging_sizes = list()
    center_sizes = list()

    for i_left, i, i_right in leading_window(window, num_genes):
        if i_right >= num_genes:        # Avoid IndexError
            i_right -= num_genes

        right_coord = stops[i_right]
        left_coord = starts[i_left]

        if right_coord < left_coord:    # Avoid negative nucleotide distance
            right_coord += length

        nucl_dist = right_coord - left_coord + 1
        leading_sizes.append(round(float(nucl_dist)/window, 2))

    for i_left, i, i_right in lagging_window(window, num_genes):
        right_coord = stops[i_right]
        left_coord = starts[i_left]

        if right_coord < left_coord:
            right_coord += length

        nucl_dist = right_coord - left_coord + 1
        lagging_sizes.append(round(float(nucl_dist)/window, 2))

    for i_left, i, i_right in center_window(window, num_genes):
        if i_right >= num_genes:        # Avoid IndexError
            i_right -= num_genes

        right_coord = stops[i_right]
        left_coord = starts[i_left]

        if right_coord < left_coord:
            right_coord += length

        nucl_dist = right_coord - left_coord + 1
        center_sizes.append(round(float(nucl_dist)/window, 2))

    return leading_sizes, center_sizes, lagging_sizes


def average_strand_changes(strands, window=WINDOW):
    """
    Calculates and returns the average number of strand changes per
    gene in leading/centered/lagging windows of genes.

    :param strands: gene orientations
    :type strands: list of int
    :param window: how many genes to average over
    :type window: int
    :return: leading_changes, centered_changes, lagging_changes
    """
    num_genes = len(strands)

    leading_changes = list()
    lagging_changes = list()
    center_changes = list()

    for i_left, i, i_right in leading_window(window, num_genes):
        cursor_strand = strands[i_left]
        strand_changes = 0
        for x in range(i_left, i_right + 1):
            if x >= num_genes:
                x -= num_genes
            if strands[x] != cursor_strand:
                strand_changes += 1
                cursor_strand = strands[x]
        leading_changes.append(round(float(strand_changes)/window, 3))

    for i_left, i, i_right in lagging_window(window, num_genes):
        cursor_strand = strands[i_left]
        strand_changes = 0
        for x in range(i_left, i_right + 1):
            if strands[x] != cursor_strand:
                strand_changes += 1
                cursor_strand = strands[x]
        lagging_changes.append(round(float(strand_changes)/window, 3))

    for i_left, i, i_right in center_window(window, num_genes):
        cursor_strand = strands[i_left]
        strand_changes = 0
        for x in range(i_left, i_right + 1):
            if x >= num_genes:
                x -= num_genes
            if strands[x] != cursor_strand:
                strand_changes += 1
                cursor_strand = strands[x]
        center_changes.append(round(float(strand_changes)/window, 3))

    return leading_changes, center_changes, lagging_changes


def build_contig_dataframe(contig):
    """
    Walks the features of the given contig and calculates the feature
    data to be used for model training and/or model prediction.

    :param contig: a contig to analyze features from
    :type contig: Bio.SeqRecord.SeqRecord
    :return: feature_dict
    """
    cds_features = [ftr for ftr in contig.features if ftr.type == "CDS"]

    temp_dict = dict()

    # These are the basis for gene density and strand bias calculations
    lefts, rights, strands = list(), list(), list()
    for feature in cds_features:
        lefts.append(feature.location.start)
        rights.append(feature.location.end)
        strands.append(feature.location.strand)

    lead_len, ctr_len, lag_len = average_gene_size(lefts, rights, len(contig))
    lead_str, ctr_str, lag_str = average_strand_changes(strands)

    temp_dict["contig_id"] = [contig.id] * len(lefts)
    temp_dict["start"] = lefts
    temp_dict["stop"] = rights
    temp_dict["strand"] = strands
    temp_dict["lag_size"] = lag_len
    temp_dict["ctr_size"] = ctr_len
    temp_dict["lead_size"] = lead_len
    temp_dict["lag_strand"] = lag_str
    temp_dict["ctr_strand"] = ctr_str
    temp_dict["lead_strand"] = lead_str

    return pd.DataFrame(temp_dict)


def smooth_by_averaging(values, window_size=25):
    """
    Smooths values by averaging them with their `window_size`
    upstream and downstream neighbors.

    :param values: the values to be smoothed
    :type values: list of float
    :param window_size: number of genes up/down-stream to average over
    :type window_size: int
    :return: smoothed_values
    """
    # Adjust window_size to avoid IndexError on short lists
    window_size = min((window_size, len(values)))
    smoothed_values = list()

    for i in range(len(values)):
        local_values = list()
        min_i, max_i = i - window_size, i + window_size
        for j in range(min_i, max_i + 1):
            if j < 0:
                j += len(values)
            elif j > len(values) - 1:
                j -= len(values)
            local_values.append(values[j])
        smoothed_values.append(average(local_values))

    return smoothed_values


def predict_prophage_genes(contig, model_path=MODEL_PATH, alpha=0.25,
                           mask=None):
    """
    Calculates the gene attributes used by the model to predict
    prophage vs bacterial genes. Then uses the classifier from
    `model_path` to make those predictions.

    :param contig: the contig to make prophage predictions in
    :type contig: Bio.SeqRecord.SeqRecord
    :param model_path: path to a binary file with sklearn model inside
    :type model_path: pathlib.Path
    :param alpha: probability above which to keep prophage prediction
    :type alpha: float
    :param mask: bitwise and will mask known/theorized bacterial genes
    :type mask: list of int
    :return: predictions
    """
    dataframe = build_contig_dataframe(contig.record)

    lead_df = pd.DataFrame()
    lead_df["ctr_size"] = dataframe.loc[:, "lead_size"]
    lead_df["ctr_strand"] = dataframe.loc[:, "lead_strand"]

    center_df = dataframe.loc[:, ["ctr_size", "ctr_strand"]]

    lag_df = pd.DataFrame()
    lag_df["ctr_size"] = dataframe.loc[:, "lag_size"]
    lag_df["ctr_strand"] = dataframe.loc[:, "lag_strand"]

    model_reader = open(model_path, "rb")
    classifier = pickle.load(model_reader)
    model_reader.close()

    lead_p = classifier.predict(lead_df)
    center_p = classifier.predict(center_df)
    lag_p = classifier.predict(lag_df)

    predictions = list()
    for x, y, z in zip(lead_p, center_p, lag_p):
        predictions.append(max((x, y, z)))

    # Store model predictions within the contig option
    contig.update_model_scores(predictions)

    predictions = smooth_by_averaging(predictions, window_size=10)

    if mask:
        for gene_i in range(len(predictions)):
            predictions[gene_i] = predictions[gene_i] * mask[gene_i]

    predictions = smooth_by_averaging(predictions, window_size=5)

    return [x >= alpha for x in predictions]


def predict_prophage_coords(contig, extend_by=0, mask=None):
    """
    Predicts prophage genes on the contig, then tries to approximate
    the coordinates associated with phage <-> bacterial transitions.

    :param contig: the contig to scan for possible prophages
    :type contig: Bio.SeqRecord.SeqRecord
    :param extend_by: number of basepairs to overextend prophages by
    :type extend_by: int
    :param mask: bitwise and will mask known/theorized bacterial genes
    :type mask: list of int
    :return: prophage_coords
    """
    gene_predictions = predict_prophage_genes(contig, mask=mask)

    prophage_coords = list()
    left, right = None, None
    n_inflections = 0
    previous_state = gene_predictions[-1]

    for i, current_state in enumerate(gene_predictions):
        if previous_state == BACTERIA and current_state == PROPHAGE:
            left = contig.genes[i].location.start - extend_by
            n_inflections += 1
        elif previous_state == PROPHAGE and current_state == BACTERIA:
            right = contig.genes[i - 1].location.end + extend_by
            prophage_coords.append((left, right))
            left, right = None, None
            n_inflections += 1
        previous_state = current_state
    # Record the last one found
    if left or right:
        prophage_coords.append((left, right))

    # Deal with edge cases (3 major possibilities):
    # Case 0: whole contig is a prophage - last gene processed was a predicted
    # prophage gene, and we had no PROPHAGE <-> BACTERIAL inflections
    if previous_state == PROPHAGE and n_inflections == 0:
        left, right = 0, len(contig.seq)
        prophage_coords.append((left, right))

    # Case 1: there's a prophage that falls off the left edge of a contig -
    # len(prophage_coords) >= 1, and no left coord for first coordinate tuple
    if len(prophage_coords) >= 1 and not prophage_coords[0][0]:
        left, right = 0, prophage_coords[0][-1]
        prophage_coords[0] = (left, right)

    # Case 2: there's a prophage that falls off the right edge of a contig -
    # len(prophage_coords) >= 1, and no right coord for last coordinate tuple
    if len(prophage_coords) >= 1 and not prophage_coords[-1][-1]:
        left, right = prophage_coords[-1][0], len(contig.seq) - 1
        prophage_coords[-1] = (left, right)

    return prophage_coords
