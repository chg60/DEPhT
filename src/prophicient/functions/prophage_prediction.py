import pickle

import pandas as pd

from prophicient import PACKAGE_DIR
from prophicient.functions.statistics import average

MODEL_PATH = PACKAGE_DIR.joinpath("data/prophage_model.pickle")

BACTERIA = 0        # Gene prediction state - bacterial
PROPHAGE = 1        # Gene prediction state - prophage


def calculate_gene_density(lefts, rights, length, windows=(5, 10, 25)):
    """
    Calculates genes per kilobase over sliding windows of different
    numbers of genes to get a robust estimate of local gene density
    for each gene.

    :param lefts: gene left coordinates
    :type lefts: list of int
    :param rights: gene right coordinates
    :type rights: list of int
    :param length: contig length
    :type length: int
    :param windows: windows to calculate gene density across
    :type windows: tuple of int
    :return: gene_densities
    """
    num_genes = len(lefts)
    gene_densities = dict()
    for window_size in windows:
        gene_densities[window_size] = list()

    for i in range(num_genes):
        for window_size in windows:
            i_left, i_right = i - window_size, i + window_size

            # Make sure indices are in range(0, num_genes)
            if i_left < 0:
                i_left += num_genes
            if i_right > num_genes - 1:
                i_right -= num_genes

            # Get left gene's left coordinate, right gene's right coordinate
            left, right = lefts[i_left], rights[i_right]
            if right < left:
                right += length

            numerator = right - left + 1
            denominator = 2 * window_size + 1

            gene_densities[window_size].append(float(numerator)/denominator)

    return gene_densities


def calculate_strand_agreement(orientations, windows=(5, 10, 25)):
    """
    Calculates percent strand agreement over sliding windows of
    different numbers of genes to get a robust estimate of local
    strand bias for each gene

    :param orientations: per gene orientations (+1/-1)
    :type orientations: list
    :param windows: window sizes to calculate pct strand agreement over
    :type windows: tuple of int
    :return: strand_agreements
    """
    num_genes = len(orientations)
    strand_agreements = dict()
    for window_size in windows:
        strand_agreements[window_size] = list()

    for i, cursor_orient in enumerate(orientations):
        for window_size in windows:
            window_orients = list()

            i_left, i_right = i - window_size, i + window_size
            for index in range(i_left, i_right + 1):
                if index < 0:
                    index += num_genes
                elif index > num_genes - 1:
                    index -= num_genes
                window_orients.append(orientations[index])

            numerator = window_orients.count(cursor_orient)
            denominator = 2 * window_size + 1

            strand_agreements[window_size].append(float(numerator)/denominator)

    return strand_agreements


def calculate_feature_dict(contig):
    """
    Walks the features of the given contig and calculates the feature
    data to be used for model training and/or model prediction.

    :param contig: a contig to analyze features from
    :type contig: Bio.SeqRecord.SeqRecord
    :return: dataframe
    """
    cds_features = [ftr for ftr in contig.features if ftr.type == "CDS"]

    gene_dict = {}

    # These are the basis for gene density and strand bias calculations
    lefts, rights, strands = list(), list(), list()
    for j, feature in enumerate(cds_features):
        left, right = feature.location.start, feature.location.end
        strand = feature.location.strand
        lefts.append(left), rights.append(right), strands.append(strand)

    # Calculate gene density attributes (in rolling windows)
    gene_densities = calculate_gene_density(lefts, rights, len(contig))
    for key, values in gene_densities.items():
        gene_dict[f"density{key}"] = values

    # Calculate strand bias attributes (in rolling windows)
    strand_biases = calculate_strand_agreement(strands)
    for key, values in strand_biases.items():
        gene_dict[f"strand{key}"] = values

    return gene_dict


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
                           beta=0.05, mask=None, iterations=3):
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
    :return: predictions
    """
    feature_dict = calculate_feature_dict(contig)
    dataframe = pd.DataFrame(feature_dict)

    with model_path.open("rb") as model_reader:
        classifier = pickle.load(model_reader)

    predictions = [x[1] for x in classifier.predict_proba(dataframe)]
    
    for _ in range(iterations):
        predictions = smooth_by_averaging(predictions, window_size=25)

        if mask is not None:
            # Amplify gene predictions with conserved bacterial genes
            for gene_i in range(len(predictions)):
                predictions[gene_i] = (predictions[gene_i] * mask[gene_i])

    if mask is not None:
        predictions = smooth_by_averaging(predictions, window_size=3)

        return [x >= beta for x in predictions]
    else:
        return [x >= alpha for x in predictions]


def predict_prophage_coords(contig, extend_by=0, mask=None):
    """
    Predicts prophage genes on the contig, then tries to approximate
    the coordinates associated with phage <-> bacterial transitions.

    :param contig: the contig to scan for possible prophages
    :type contig: Bio.SeqRecord.SeqRecord
    :param extend_by: number of basepairs to overextend prophages by
    :type extend_by: int
    :return: prophage_coords
    """
    gene_predictions = predict_prophage_genes(contig, mask=mask) 

    prophage_coords = list()
    left, right = None, None
    n_inflections = 0
    previous_state = gene_predictions[-1]

    # predictions are only made on the cds features
    cds_features = [feat for feat in contig.features if feat.type == "CDS"]

    for i, current_state in enumerate(gene_predictions):
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
