from depht.functions.prophage_prediction import (build_contig_dataframe,
                                                 smooth_by_averaging,
                                                 filter_prophage_signal,
                                                 predict_prophage_genes,
                                                 BACTERIA, PROPHAGE)


def get_genomic_islands(contig, classifier, extend_by=0, mask=None):
    gene_predictions = predict_genomic_island_genes(contig, classifier, mask)

    prophage_coords = list()
    left, right = None, None
    n_inflections = 0
    previous_state = gene_predictions[0]

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


def predict_genomic_island_genes(contig, classifier, predictions, alpha=0.25):
    predictions = smooth_by_averaging(predictions, window_size=5)
    phage_gene_mask = predict_prophage_genes(contig, classifier) 

    for i in range(len(phage_gene_mask)):
        phage_gene_mask[i] = int(not phage_gene_mask[i])

    for gene_i in range(len(predictions)):
        predictions[gene_i] = predictions[gene_i] * phage_gene_mask[gene_i]

    predictions = smooth_by_averaging(predictions, window_size=5)

    island_signal = [x >= alpha for x in predictions]

    return island_signal 

