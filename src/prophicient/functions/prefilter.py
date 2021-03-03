from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (FeatureLocation, SeqFeature)

from numpy import (average, std)
from Prophicient.functions import att
from Prophicient.utilities import plot_gene_density


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULTS = {"bin_width": plot_gene_density.BIN_WIDTH,
            "window_size": plot_gene_density.WINDOW_SIZE,
            "F": -0.5, "C": 2}


def main(bin_width=DEFAULTS["bin_width"], window_size=DEFAULTS["window_size"],
         F=DEFAULTS["F"], C=DEFAULTS["C"]):
    record = SeqIO.read("GD43A.gbk", "gb")

    gene_dense_coords = prefilter_genome(record, window_size=window_size,
                                         F=F, C=C)

    gene_dense_regions = list()
    for coords in gene_dense_coords:
        feature = SeqFeature(FeatureLocation(coords[0], coords[1]),
                             strand=1)

        gene_dense_regions.append(feature)

    phage_regions = list()
    for region in gene_dense_regions:
        sequence = str(region.extract(record.seq))

        half = len(sequence) // 2
        l_region = sequence[:half-1]
        r_region = sequence[half+1:]

        print("Found elevated gene dense region "
              f"{len(sequence) // 1000} kb long: "
              f"({region.location.start}, {region.location.end})")

        attL_feature, attR_feature = att.find_attatchment_site(
                                                        l_region, r_region)
        att_seq = str(attL_feature.extract(Seq(l_region)))

        prophage_start = region.location.start + attL_feature.location.start
        prophage_end = (region.location.start + half + 1 +
                        attR_feature.location.end)
        phage_regions.append((prophage_start, prophage_end))

        print("\tPutative prophage "
              f"{(prophage_end - prophage_start) // 1000} kb long @ "
              f"({prophage_start}, {prophage_end}) "
              f"with att: {att_seq}")

    return phage_regions


def prefilter_genome(record, bin_width=DEFAULTS["bin_width"],
                     window_size=DEFAULTS["window_size"],
                     F=DEFAULTS["F"], C=DEFAULTS["C"]):
    nu_map = plot_gene_density.count_genes_per_interval(
                                                record, interval=bin_width)

    rho_map = calculate_rho_map(len(record), nu_map, window_size,
                                bin_width)

    rho_data = [rho for rho in rho_map.values()]

    rho_avg = average(rho_data)
    rho_std = std(rho_data)
    floor = rho_avg - (rho_std * F)

    contig_rho = compile_rho_contigs(rho_map, floor)

    ceiling = rho_avg + (rho_std * C)

    target_contig_rho = list()
    for contig in contig_rho:
        if contig[0][1] >= ceiling:
            contig.sort(key=lambda x: x[0])
            target_contig_rho.append(contig)

    gene_dense_coords = list()
    for contig in target_contig_rho:
        start = contig[0][0] - window_size // 2
        end = contig[-1][0] + window_size // 2

        gene_dense_coords.append((start, end))

    return gene_dense_coords


def calculate_rho_map(record_len, nu_map, window_size, bin_width):
    rho_dict = dict()
    for x in range(0, record_len, bin_width):
        gene_count = 0
        for j in range(x - window_size // 2, x + window_size // 2, bin_width):
            gene_count += nu_map[j]

        rho = gene_count / window_size * bin_width
        rho_dict[x] = rho

    return rho_dict


def compile_rho_contigs(rho_dict, floor):
    rho_indicies = [rho_tuple for rho_tuple in rho_dict.items()]

    rho_indicies.sort(key=lambda x: x[0])
    num_rho = len(rho_indicies)

    index_rho_map = dict()
    index_position_map = dict()
    for i in range(len(rho_indicies)):
        index_position_map[i] = rho_indicies[i][0]

        rho = rho_indicies[i][1]
        if rho >= floor:
            index_rho_map[i] = rho

    contig_rho = list()
    contig = None
    for i in range(num_rho):
        rho = index_rho_map.get(i)

        if rho is not None:
            if contig is None:
                contig = list()

            contig.append((index_position_map[i], rho))

            continue

        if contig is not None:
            contig.sort(key=lambda x: x[1], reverse=True)
            contig_rho.append(contig)

        contig = None

    return contig_rho


if __name__ == "__main__":
    main()
