# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULTS = {"bin_width": 1000, "window_size": 50000}


def count_genes_per_interval(record, interval=DEFAULTS["bin_width"]):
    """
    Function to count the number of genes per `interval` nucleotides in the
    record.

    Counts gene starts, but could easily be adapted to work in some other way

    :param record: a Genbank record to process
    :type record: Bio.SeqRecord.SeqRecord
    :param interval: the number of bases to advance when indexing
    :type interval: int
    :return: genes_per_interval
    """

    genes_per_interval = dict()
    max_interval = len(record) // interval * interval

    # Initialize every interval's counter to 0
    for interval_index in range(0, len(record), interval):
        genes_per_interval[interval_index] = 0

    # Iterate over genes, and place them in their appropriate bin
    for gene in record.features:
        interval_index = gene.location.start // interval * interval
        genes_per_interval[interval_index] += 1

    # Wrap around ends of the genome
    for interval_index in range(0, DEFAULTS["window_size"], interval):
        genes_per_interval[interval_index + max_interval] = genes_per_interval[interval_index]
    for interval_index in range(max_interval - DEFAULTS["window_size"],
                                max_interval, interval):
        genes_per_interval[interval_index - max_interval] = genes_per_interval[interval_index]

    return genes_per_interval


def calculate_density_map(record_len, gene_count_map, window_size, bin_width):
    index_rho_dict = dict()
    for x in range(0, record_len, bin_width):
        gene_count = 0
        for j in range(x - window_size // 2, x + window_size // 2, bin_width):
            gene_count += gene_count_map[j]

        rho = gene_count / window_size * bin_width
        index_rho_dict[x] = rho

    return index_rho_dict


def compile_region_contigs(rho_dict, eps):
    rho_indicies = [rho_tuple for rho_tuple in rho_dict.items()]

    rho_indicies.sort(key=lambda x: x[0])
    num_rho = len(rho_indicies)

    index_rho_map = dict()
    index_position_map = dict()
    for i in range(len(rho_indicies)):
        index_position_map[i] = rho_indicies[i][0]

        rho = rho_indicies[i][1]
        if rho >= eps:
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
