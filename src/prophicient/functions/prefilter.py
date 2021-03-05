from Bio.SeqFeature import (FeatureLocation, SeqFeature)

from numpy import (average, std)
from Prophicient.functions import (att, gene_density)


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULTS = {"bin_width": gene_density.DEFAULTS["bin_size"],
            "window_size": gene_density.DEFAULTS["window_size"],
            "F": -0.5, "C": 2}


def prefilter_genome(record, bin_width=DEFAULTS["bin_width"],
                     window_size=DEFAULTS["window_size"], F=DEFAULTS["F"],
                     C=DEFAULTS["C"]):
    gene_dense_regions = get_gene_dense_regions(
                                            record, window_size=window_size,
                                            F=F, C=C)

    phage_regions = list()
    for region in gene_dense_regions:
        sequence = str(region.extract(record.seq))

        half = len(sequence) // 2
        l_region = sequence[:half-1]
        r_region = sequence[half+1:]

        attL_feature, attR_feature = att.find_attatchment_site(
                                                        l_region, r_region)

        prophage_start = region.location.start + attL_feature.location.start
        prophage_end = (region.location.start + half + 1 +
                        attR_feature.location.end)

        phage_regions.append((prophage_start, prophage_end))

    return phage_regions


def get_gene_dense_regions(record, bin_width=DEFAULTS["bin_width"],
                           window_size=DEFAULTS["window_size"],
                           F=DEFAULTS["F"], C=DEFAULTS["C"]):
    pos_nu_map = gene_density.count_genes_per_interval(
                                                    record, interval=bin_width)

    index_rho_map = gene_density.calculate_density_map(
                                        len(record), pos_nu_map, window_size,
                                        bin_width)

    rho_data = [rho for rho in index_rho_map.values()]

    rho_avg = average(rho_data)
    rho_std = std(rho_data)
    eps = rho_avg - (rho_std * F)

    region_contigs = gene_density.compile_region_contigs(index_rho_map, eps)

    ceiling = rho_avg + (rho_std * C)

    target_contig_rho = list()
    for contig in region_contigs:
        if contig[0][1] >= ceiling:
            contig.sort(key=lambda x: x[0])
            target_contig_rho.append(contig)

    gene_dense_coords = list()
    for contig in target_contig_rho:
        start = contig[0][0] - window_size // 2
        end = contig[-1][0] + window_size // 2

        gene_dense_coords.append((start, end))

    gene_dense_regions = list()
    for coords in gene_dense_coords:
        feature = SeqFeature(FeatureLocation(coords[0], coords[1]),
                             strand=1)

        gene_dense_regions.append(feature)

    return gene_dense_regions
