import math

from networkx import DiGraph
from scipy.stats import zscore

from depht.classes import kmers
from depht.classes.prophage import DEFAULT_PRODUCT
from depht.functions.blastn import blastn
from depht.functions.fasta import write_fasta
from depht.functions.statistics import transform

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
KMER_SIZE = 5
MIN_ATT_SCORE = 2.3
EVALUE_FILTER = 10000

L_SEQ_NAME = "putative_attL_region"
R_SEQ_NAME = "putative_attR_region"

AQ_WEIGHT = 1
IP_WEIGHT = 0.6
MC_WEIGHT = 0.9
TR_WEIGHT = 0
RC_WEIGHT = 1

DEFAULTS = {"k": 5, "fpp": 0.0001, "outfmt": 10}


def find_attachment_site(prophage, l_seq, r_seq,
                         reference_db_path, tmp_dir, sort_key,
                         k=KMER_SIZE, min_score=MIN_ATT_SCORE,
                         l_name=L_SEQ_NAME, r_name=R_SEQ_NAME):
    """Given the sequences of a putative attL region and putative attR region,
    find the most probable attachment site, dictated by the sequence's length
    and it's distance from the predicted origin position.

    :param prophage: Prophage object to find an attachment site for
    :type prophage: prophicient.classes.prophage.Prophage
    :param l_seq: The sequence of a putative attL region.
    :type l_seq: str
    :param r_seq: The sequence of a putative attR region.
    :type r_seq: str
    :param tmp_dir: The directory for which to write sequences and files to
    :type tmp_dir: pathlib.Path
    :param k: Length of the word size storerd in the DeBruijn graph.
    :type k: int
    :param l_name: Name to give to the putative attL region sequence.
    :type l_name: str
    :param r_name: Name to give to the putative attR region sequence.
    :type r_name: str
    :return: A tuple of information associated with the detected att site.
    :rtype: tuple
    """
    # Write the putative attL region to file
    l_seq_path = tmp_dir.joinpath(f"{l_name}.fasta")
    write_fasta([l_name], [l_seq], l_seq_path)

    # Write the putative attR region to file
    r_seq_path = tmp_dir.joinpath(f"{r_name}.fasta")
    write_fasta([r_name], [r_seq], r_seq_path)
    # Calculate and store right region's coordinate start for easy access
    r_seq_start = prophage.end - len(r_seq)

    paired_ref_map = find_reference_att_sites(l_seq_path, r_seq_path,
                                              reference_db_path, tmp_dir,
                                              k, sort_key,
                                              prophage.start, r_seq_start)

    # Use BLASTn to retrieve putative attachment site sequences
    kmer_contigs = blast_attachment_site(l_seq_path, r_seq_path, tmp_dir, k=k,
                                         evalue=len(l_seq))

    if not kmer_contigs:
        return

    transform_kmer_contig_bitscores(kmer_contigs)

    # Score putative attachment site sequences
    # TODO paramaterize and allow access to change the 'exponent' variable
    scored_kmer_contigs = []
    for kmer_contig in kmer_contigs:
        scores = score_kmer(kmer_contig, prophage, paired_ref_map, r_seq_start)

        scored_kmer_contigs.append((kmer_contig, scores))

    # Sort attachment site sequences by score
    scored_kmer_contigs.sort(key=lambda x: x[1][0], reverse=True)

    att_table_path = tmp_dir.joinpath("att.txt")
    dump_attachment_sites(prophage, scored_kmer_contigs, att_table_path,
                          r_seq_start)

    kmer_contig, scores = scored_kmer_contigs[0]

    if scores[0] < min_score:
        return

    new_start = prophage.start + kmer_contig[1]
    new_end = r_seq_start + kmer_contig[2]

    return new_start, new_end, scores, kmer_contig[0]


# MAIN HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def find_reference_att_sites(left_seq_path, right_seq_path, reference_db_path,
                             tmp_dir, k, sort_key, l_seq_start, r_seq_start):
    # BLASTn both regions against the reference database
    left_map = build_reference_map(left_seq_path, reference_db_path, tmp_dir)
    right_map = build_reference_map(right_seq_path, reference_db_path, tmp_dir)

    ref_ids = list(set(left_map.keys()).intersection(set(right_map.keys())))

    paired_ref_map = pair_reference_maps(ref_ids, left_map, right_map,
                                         k, sort_key, l_seq_start, r_seq_start)

    return paired_ref_map


def blast_attachment_site(l_seq_path, r_seq_path, tmp_dir, k=KMER_SIZE,
                          evalue=EVALUE_FILTER):
    """Given the path to files containink the putative attL region and
    putative attR region, BLASTn the sequences of both regions against each
    other and retrieve matching sequences and their positions.

    :param l_seq_path: path to the sequence of the putative attL region
    :type l_seq_path: pathlib.Path
    :param r_seq_path: path to the sequence of the putative attR region
    :type r_seq_path: pathlib.Path
    :param tmp_dir: path where temporary files can go
    :type tmp_dir: pathlib.Path
    :param k: Length of the word size used by the BLAST algorithm
    :type k: int
    :return: A list of contigs and their positions in the sequence and graph.
    :rtype: list(tuple(str, int, int))
    """
    blast_results = blastn(
        l_seq_path, r_seq_path, tmp_dir, mode="subject", word_size=k,
        evalue=evalue, gapopen=10, gapextend=4)

    kmer_contigs = []
    for result in blast_results:
        # Append the contig and its positions to the list
        kmer_contig = [result["qseq"],
                       int(result["qstart"]) - 1, int(result["send"]),
                       float(result["bitscore"])]
        kmer_contigs.append(kmer_contig)

    return kmer_contigs


def graph_attachment_site(l_seq, r_seq, k=KMER_SIZE):
    """Given the sequences of a putative attL region and putative attR region,
    kmer count the sequences of both regions using a DeBruijn graph that has
    been optimized by reducing the inputted words with a Bloom Filter to
    find sequence contigs and their positions.

    :param l_seq: The sequence of a putative attL region.
    :type l_seq: str
    :param r_seq: The sequence of a putative attR region.
    :type r_seq: str
    :param k: Length of the word size stored in the DeBruijn graph.
    :type k: int
    :return: A list of contigs and their positions in the sequence and graph.
    :rtype: list(tuple(str, int, int))
    """
    # Load a Bloom Filter from the sequence of the putative attL region
    bfilter = load_bfilter(l_seq, k=k)
    # Initialize a DeBruijn graph from the sequence of the putative attR region
    # limited by the kmers represented in the loaded Bloom Filter
    deb_graph = create_debruijn_graph(bfilter, r_seq)
    # Retrieve kmer contigs by tracing the DeBruijn graph and
    # the sequence of the putative attL region
    kmer_contigs = traverse_debruijn_graph(deb_graph, l_seq)

    return kmer_contigs


def dump_attachment_sites(prophage, scored_kmer_contigs, outpath, r_seq_start):
    filehandle = outpath.open(mode="w")

    for kmer_contig, scores in scored_kmer_contigs:
        new_start = prophage.start + kmer_contig[1]
        new_end = r_seq_start + kmer_contig[2]

        att_line_data = [new_start, new_end, len(kmer_contig[0])]
        score_line_data = [round(score, 2) for score in scores]
        seq_data = [kmer_contig[0]]

        line_data = att_line_data + score_line_data + seq_data
        line_data = [str(line_entry) for line_entry in line_data]

        filehandle.write("\t".join(line_data))
        filehandle.write("\n")

    filehandle.close()


# REFERENCE BLASTING HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def build_reference_map(sequence_path, reference_db_path, tmp_dir):
    """Maps sequence BLASTn aligned reference genome IDs to their respective
    alignment result data.

    :param sequence_path: Path to query to be aligned to the reference database
    :type sequence_path: pathlib.Path
    :param reference_db_path: Path to the database of references to search
    :type reference_db_path: pathlib.Path
    :param tmp_dir: Working directory to place BLASTn inputs and outputs
    :type tmp_dir: pathlib.Path
    :return: A map of aligned reference genome IDs to alignment result data
    """
    # Try to retrieve reference results for the sequence to the references
    blast_results = blastn(sequence_path, reference_db_path, tmp_dir)

    reference_map = dict()
    for blast_result in blast_results:
        # Checks to see if the sequence reference ID has already been stored
        results = reference_map.get(blast_result["sseqid"], list())

        results.append(blast_result)
        reference_map[blast_result["sseqid"]] = results

    return reference_map


def pair_reference_maps(ref_ids, left_map, right_map, k, sort_key,
                        l_seq_start, r_seq_start):
    paired_ref_map = {}
    for ref_id in ref_ids:
        ref_data = []

        for l_data in left_map[ref_id]:
            for r_data in right_map[ref_id]:
                # Find the coordinate ranges of the aligned reference genomes
                left_ref_range = range(int(l_data["sstart"]),
                                       int(l_data["send"]))
                right_ref_range = range(int(r_data["sstart"]),
                                        int(r_data["send"]))

                # Determine if the coordinate ranges overlap
                overlap_range = set(left_ref_range).intersection(
                                                        set(right_ref_range))
                overlap_len = len(overlap_range)

                # If the overlap length meets the minimum att length
                # treat the sequence as a putative attB that we can use
                # as a reference to set the boundaries for the prophage
                if overlap_len >= k:
                    # Find the right coordinate of the putative attL
                    # in the aligned left region of the prophage
                    l_qend = int(l_data["qend"])

                    # Find the left coordinate of the putative attR
                    # in the aligned right region of the prophage
                    r_qstart = int(r_data["qstart"])

                    new_start = (l_seq_start + int(l_qend) - overlap_len)
                    new_end = (r_seq_start + int(r_qstart) + overlap_len)

                    score = float(l_data[sort_key]) + float(r_data[sort_key])
                    att_data = (new_start, new_end, overlap_len, score)

                    ref_data.append(att_data)

            ref_data.sort(key=lambda x: x[3], reverse=True)

            if ref_data:
                paired_ref_map[ref_id] = ref_data[0]

    return paired_ref_map


# SCORING HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def transform_kmer_contig_bitscores(kmer_contigs):
    """Transform the kmer contig bitscores into z-scores
    :param kmer_contigs: A list of kmers and their positions in the sequence
    :type kmer_contigs: list
    """
    bitscores = []
    for kmer_contig in kmer_contigs:
        bitscores.append(kmer_contig[3])

    zscores = zscore(bitscores)
    transform(zscores, min_t=0, max_t=1)

    for i, kmer_contig in enumerate(kmer_contigs):
        kmer_contig.append(zscores[i])


def score_kmer(kmer_contig, prophage, paired_ref_map, r_seq_start):
    """Score kmer contigs with a complete holsitic approach.

    :param kmer_contig: A tuple containing the contig, and its positions
    :type kmer_contig: tuple(str, int, int)
    :return: The score of the given kmer
    :rtype: float
    """
    attL_pos = prophage.start + kmer_contig[1]
    attR_pos = r_seq_start + kmer_contig[2]

    att_quality_score = score_att_quality(kmer_contig[4])

    int_proximity_score, int_distance = score_integrase_proximity(
                                                    prophage,
                                                    attL_pos - prophage.start,
                                                    attR_pos - prophage.start)

    model_cov_score, model_coverage = score_model_coverage(attR_pos - attL_pos,
                                                           len(prophage.seq))

    reference_score, reference_bitscore = score_reference_concurrence(
                                        attL_pos, attR_pos,
                                        len(kmer_contig[0]), paired_ref_map)

    composite_score = ((att_quality_score + int_proximity_score +
                        model_cov_score) *
                       (3 / (AQ_WEIGHT + IP_WEIGHT + MC_WEIGHT)) +
                       reference_score)

    return (composite_score, att_quality_score, kmer_contig[3],
            int_proximity_score, int_distance,
            model_cov_score, model_coverage,
            reference_score, reference_bitscore)


def score_att_quality(normalized_bitscore, weight=AQ_WEIGHT):
    score = normalized_bitscore

    if score < 0:
        score = 0.0
    elif score > 1:
        score = 1.0

    weighted_score = score * weight
    return weighted_score


def score_integrase_proximity(prophage, attL_pos, attR_pos, base_dist=1500,
                              weight=IP_WEIGHT):
    int_dist = None
    for feature in prophage.record.features:
        if feature.type != "CDS":
            continue

        product = feature.qualifiers.get("product", [DEFAULT_PRODUCT])[0]

        if "integrase" in product:
            left_int_dist = int(feature.location.start - attL_pos)
            right_int_dist = int(attR_pos - feature.location.end)

            if int_dist is None:
                int_dist = left_int_dist
            elif int_dist < 0 and left_int_dist > int_dist:
                int_dist = left_int_dist
            else:
                if left_int_dist < int_dist and left_int_dist > 0:
                    int_dist = left_int_dist

            if int_dist is None:
                int_dist = right_int_dist
            elif int_dist < 0 and right_int_dist > int_dist:
                int_dist = right_int_dist
            else:
                if right_int_dist < int_dist and right_int_dist > 0:
                    int_dist = right_int_dist

    if int_dist is None:
        score = -1
        int_dist = -1
    elif int_dist < 0:
        score = -1
    elif int_dist == 0:
        score = 1
    else:
        score = (base_dist / int_dist) ** 2

        if score > 1:
            score = 1

    weighted_score = score * weight
    return weighted_score, int_dist


def score_trna_overlap(prophage, attL_pos, attR_pos, att_len,
                       weight=TR_WEIGHT):
    overlap = 0
    for feature in prophage.record.features:
        if feature.type != "tRNA":
            continue

        trna_range = set(range(feature.location.start, feature.location.end))

        attL_range = set(range(attL_pos, attL_pos + att_len))
        attR_range = set(range(attR_pos - att_len, attR_pos))

        if trna_range.intersection(attL_range):
            overlap = 1
        elif trna_range.intersection(attR_range):
            overlap = 1

        if overlap > 0:
            break

    if overlap <= 0:
        return (0, overlap)
    else:
        return (overlap * TR_WEIGHT, overlap)


def score_model_coverage(putative_len, model_len, weight=MC_WEIGHT):
    score = putative_len / model_len

    if score > 1:
        score = 1.0

    weighted_score = score * weight
    return score, weighted_score


def score_reference_concurrence(attL_pos, attR_pos, att_len, paired_ref_map,
                                base=10, weight=RC_WEIGHT):
    """Return a score dependant on the coordinate range overlap of the
    kmer contig with a preducted reference attB.
    """
    ref_bitscore = None
    for ref_id, ref_data in paired_ref_map.items():
        attL_range = set(range(attL_pos, (attL_pos + att_len)))
        ref_attL_range = set(range(ref_data[0], (ref_data[0] + att_len)))

        if not attL_range.intersection(ref_attL_range):
            continue

        attR_range = set(range((attR_pos - att_len), attR_pos))
        ref_attR_range = set(range((ref_data[1] - att_len), ref_data[1]))

        if not attR_range.intersection(ref_attR_range):
            continue

        if ref_bitscore is None:
            ref_bitscore = ref_data[3]
        else:
            if ref_data[3] > ref_bitscore:
                ref_bitscore = ref_data[3]

    if ref_bitscore is None:
        return 0.0, 0

    if ref_bitscore < 10:
        ref_bitscore = 10

    score = 1 - (1 / math.log(ref_bitscore, base))

    if score < 0:
        score = 0.0
    elif score > 1:
        score = 1.0

    weighted_score = score * weight
    return weighted_score, ref_bitscore


# DEB-GRAPH KMER-COUNTING HELPER FUNCTIONS
# -----------------------------------------------------------------------------


def load_bfilter(sequence, k=KMER_SIZE):
    """Loads a Bloom Filter with the kmers from a given sequence.

    :param sequence: A sequence to load a Bloom Filter with.
    :type sequence: str
    :param k: Length of the word size stored in the DeBruijn graph.
    :type k: int
    :return: A Bloom Filter loaded with the kmers from the given sequence
    :rtype: prophicient.classes.kmer.BloomFilter
    """
    bfilter = kmers.BloomFilter(len(sequence))

    for pos, kmer in count_kmers(sequence, k):
        bfilter.add(kmer)

    return bfilter


def create_debruijn_graph(bfilter, sequence, k=KMER_SIZE):
    """Create a DeBruijn graph from a sequence with a loaded Bloom Filter
    as a guide.

    :param bfilter: A Bloom Filter that has been loaded with a sequence.
    :type bfilter: prophicient.classes.kmer.BloomFilter
    :param sequence: A sequence to construct a DeBruijn graph from.
    :type sequence: str
    :param k: Length of the word size to store in the DeBruijn graph.
    :type k: int
    :return: A DeBruijn graph loaded from the given sequence.
    :rtype: networkx.DiGraph
    """
    # Initialize the DiGraph
    deb_graph = DiGraph()

    # Initialize anchor position and previous kmer variables
    anchor_pos = None
    prev_kmer = None
    # Iterate over the kmers and their positions in the given sequence
    for pos, kmer in count_kmers(sequence, k):
        # If the kmer has not been seen by the Bloom Filter, continue
        if not bfilter.check(kmer):
            anchor_pos = None
            prev_kmer = None
            continue

        # If the anchor position for the current stretch of kmers has not been
        # set, set the anchor position to the current kmer position
        if anchor_pos is None:
            anchor_pos = pos

        # Try to get a node associated with the current kmer
        node = deb_graph.nodes.get(kmer)
        # If no node exists, create the node associated with the kmer
        if node is None:
            # Create a set with the current anchor position
            anchors = set()
            anchors.add(anchor_pos)

            # Find the kmer's relative position to the current anchor
            rel_pos_lookup = dict()
            rel_pos = set()
            rel_pos.add(pos - anchor_pos)
            rel_pos_lookup[anchor_pos] = rel_pos

            # Add the new node to the DeBruijn graph, along with
            # the set of anchor positions and relative positions
            deb_graph.add_node(kmer, anchors=anchors,
                               rel_pos_lookup=rel_pos_lookup)
        # If the node exists, update the node's position metadata
        else:
            # Add the current anchor position to the node's metadata
            anchors = node.get("anchors")
            anchors.add(anchor_pos)

            rel_pos_lookup = node.get("rel_pos_lookup")

            # Add the kmer's relative position to the current anchor
            # to the node's metadata
            rel_pos = rel_pos_lookup.get(anchor_pos, set())
            rel_pos.add(pos - anchor_pos)
            rel_pos_lookup[anchor_pos] = rel_pos

        # Set the previous kmer variable
        if prev_kmer is None:
            prev_kmer = kmer
            continue

        # If no edge between the previous kmer and the current kmer had existed
        # create a new edge
        edge = deb_graph[prev_kmer].get(kmer)
        if edge is None:
            deb_graph.add_edge(prev_kmer, kmer)

        prev_kmer = kmer

    return deb_graph


def traverse_debruijn_graph(deb_graph, sequence, k=KMER_SIZE):
    """Iterate over a sequence and, with a DeBruijn graph as a guide,
    find kmer contigs.

    :param deb_graph: A DeBruijn graph of words in a sequence.
    :type deb_graph: networkx.DiGraph
    :param sequence: A generator yielding tuples of kmers and their positions.
    :type sequence: str
    :param k: Length of the word size stored in the DeBruijn graph.
    :type k: int
    :return: A list of contigs and their positions in the sequence and graph.
    :rtype: list(tuple(str, int, int))
    """
    kmer_contigs = []

    # Initialize the kmer generator
    kmer_gen = count_kmers(sequence, k)
    while True:
        # Try to iterate to the next kmer and position
        try:
            pos, kmer = next(kmer_gen)
        except StopIteration:
            break

        # Get the graph node associated with the current kmer
        node = deb_graph.nodes.get(kmer)
        if node is not None:
            # Try to find and stitch together a multiple kmer contig
            graph_pos, contig = stitch_kmer_contigs(deb_graph, kmer_gen,
                                                    kmer, k)
            kmer_contigs.append((contig, pos, graph_pos))

    return kmer_contigs


def stitch_kmer_contigs(deb_graph, kmer_gen, kmer, k):
    """Stitches together kmers from a contiguous sequence represented in
    a DeBruijn graph beginning with a given kmer.

    :param deb_graph: A DeBruijn graph of words in a sequence.
    :type deb_graph: networkx.DiGraph
    :param kmer_gen: A generator yielding tuples of kmers and their positions.
    :type kmer_gen: generator
    :param kmer: The kmer word to begin the current contig.
    :type kmer: str
    :param k: Length of the word size stored in the DeBruijn graph.
    :type k: int
    :return: A tuple of the position and sequence of a DeBruijn graph contig
    :rtype: tuple(int, str)
    """
    # Intializes the contig variable
    contig = kmer
    # Retrieves the node associated with the beginning kmer
    node = deb_graph.nodes[kmer]
    # Retrieves the anchor positions associated with the beginning kmer
    graph_anchor_set = node["anchors"]
    # Retrieves the positions relative to an anchor
    # associated with the beginning kmer
    graph_rel_pos_set = build_graph_rel_pos_set(node, graph_anchor_set)
    # Retrieves the loosest (rightmost) value the kmer could appear at
    graph_pos = get_graph_pos(graph_anchor_set, graph_rel_pos_set, k)

    # Iterates until the next kmer in the generator cannot yield a proper path
    # or the generator runs out of values
    prev_kmer = kmer
    while True:
        # Breaks if the generator runs out of values
        try:
            pos, kmer = next(kmer_gen)
        except StopIteration:
            break

        # Breaks if the next kmer does not exist within the DeBruijn  graph
        node = deb_graph.nodes.get(kmer)
        if node is None:
            break

        # Test if the next kmer has ever appeared next to the currrent kmer
        # If yes, return the anchor positions they appear together on
        # with the set of the previous kmers in sequence
        graph_anchor_set = test_valid_kmer_path(deb_graph, graph_anchor_set,
                                                prev_kmer, kmer)

        # If the anchor positions for the kmer pair do not align
        # with the rest of the kmers in sequence, break
        if not graph_anchor_set:
            break

        # If the relative positions for the kmer pair indicate that
        # the the kmers cannot exist adjacent to each other with the rest
        # of the kmers in sequence, break
        graph_rel_pos_set, graph_anchor_set = test_valid_kmer_position(
                                                 deb_graph, graph_rel_pos_set,
                                                 graph_anchor_set, kmer)
        if (not graph_rel_pos_set) or (not graph_anchor_set):
            break

        # Update kmer, contig, and loosest graph position
        prev_kmer = kmer
        contig += kmer[-1]
        graph_pos = get_graph_pos(graph_anchor_set, graph_rel_pos_set, k)

    return graph_pos, contig


def test_valid_kmer_path(deb_graph, graph_anchor_set, prev_kmer, kmer):
    """Determines if a kmer has ever appeared adjacent to a previous kmer
    in a DeBruijn graph.

    :param deb_graph: A DeBruijn graph of words in a sequence.
    :type deb_graph: networkx.DiGraph
    :param graph_anchor_set: Positions that the contig could be anchored at.
    :type graph_anchor_set: set(int)
    :param prev_kmer: The previous kmer in a DeBruijn graph.
    :type prev_kmer: str
    :param kmer: The kmer to append to the current contig.
    :type kmer: str
    :return: The positions that the kmer's contigs could appear at
    :rtype: set(int)
    """
    node = deb_graph.nodes[kmer]

    # Tests if the kmer has ever appeared adjacent to the previous kmer
    edge = deb_graph[prev_kmer].get(kmer)
    if edge is None:
        return

    # Retrieves the graph contig anchor positions that this edge belongs to
    graph_anchor_set = graph_anchor_set.intersection(node["anchors"])
    return graph_anchor_set


def test_valid_kmer_position(deb_graph, graph_rel_pos_set, graph_anchor_set,
                             kmer):
    """Tests whether the kmer can be appended to the given contig
    by using the relative positions that the kmer appears in the
    DeBruijn graph.

    :param deb_graph: A DeBruijn graph of words in a sequence.
    :type deb_graph: networkx.DiGraph
    :param graph_rel_pos_set: Positions that the previous kmer could appear at.
    :type graph_rel_pos_set: set(int)
    :param graph_anchor_set: Positions that the contig could be anchored at.
    :type graph_anchor_set: set(int)
    :param kmer: The kmer word to append to the current contig.
    :type kmer: str
    :return: The positions that the kmer and its contigs could appear at
    :rtype: tuple(set(int), set(int))
    """
    # Increments the previous kmer's valid positions
    inc_graph_rel_pos_set = increment_graph_rel_pos_set(graph_rel_pos_set)

    graph_rel_pos_anchor_set = set()
    int_graph_rel_pos_set = set()

    # Iterates over the known relative positions for the given kmer
    for anchor, rel_pos_set in deb_graph.nodes[kmer]["rel_pos_lookup"].items():
        # Finds positions where the kmer exists adjacent to the previous kmer
        for rel_pos in rel_pos_set:
            # If there is such a position, store this as well as
            # the graph contig that the relative position belongs to
            if rel_pos in inc_graph_rel_pos_set:
                graph_rel_pos_anchor_set.add(anchor)
                int_graph_rel_pos_set.add(rel_pos)

    # Find the graph contig positions that this particular ordering of kmers
    # could appear at, through set intersection of their positions
    graph_anchor_set = graph_anchor_set.intersection(graph_rel_pos_anchor_set)

    return int_graph_rel_pos_set, graph_anchor_set


def build_graph_rel_pos_set(node, graph_anchor_set):
    """Builds a set of all of the relative positions that a kmer appears in
    with respect to the position of the start of their contig.

    :param node: A DeBruijn graph node.
    :type node: networkx.Graph.node
    :param graph_anchor_set: Positions of the contigs that contain the kmer.
    :type graph_anchor_set: set(int)
    :return: A set of all of the relative positions a kmer appears in
    :rtype: set(int)
    """
    graph_rel_pos_set = set()

    # Iterates over all the DeBruij graph contigs that the kmer node appears in
    for graph_anchor in graph_anchor_set:
        # Finds all the relative positions where the kmer appears in the contig
        node_rel_pos_set = node["rel_pos_lookup"][graph_anchor]

        # Adds each relative position to the set
        for rel_pos in node_rel_pos_set:
            graph_rel_pos_set.add(rel_pos)

    return graph_rel_pos_set


def increment_graph_rel_pos_set(graph_rel_pos_set):
    """Takes a set of DeBruijn graph word positions and increments each
    position by one.

    :param graph_rel_pos_set: A set of DeBruijn graph word positions.
    :type graph_rel_pos_set: set(int)
    :returns: A deep copy of the positions of each word, incremented by one.
    :rtype: set(int)
    """
    inc_graph_rel_pos_set = set()

    for rel_pos in graph_rel_pos_set:
        inc_graph_rel_pos_set.add(rel_pos + 1)

    return inc_graph_rel_pos_set


def get_graph_pos(graph_anchor_set, graph_rel_pos_set, k):
    """Returns the loosest value permitted by the relative positions stored
    within a nucleotide sequence DeBruijn graph.

    :param graph_anchor_set: Position of the beginning of related graph contigs
    :type graph_anchor_set: set(int)
    :param graph_rel_pos_set: Positions of the beginning of a graph kmer (node)
    :type graph_anchor_set: set(int)
    :param k: Length of the word size stored in the DeBruijn graph
    :type k: int
    :rtype: int
    """
    # Returns the right-most value the kmer could appear at.
    return max(graph_anchor_set) + max(graph_rel_pos_set) + k


def count_kmers(sequence, k):
    """A generator function that takes a sequence and splits it into word
    sizes of k and the index of the word with respect to the sequence.

    :param sequence: A sequence to split into k sized word.
    :type sequence: str
    :param k: The size of words to split the specified sequence into.
    :type k: int
    :returns: Generator that yields k length words from the sequence.
    """
    # If the sequence is not long enough, return a faux generator result
    if k >= len(sequence):
        return [(0, sequence)]

    # Generator that returns k-mers and their 0-indexed position
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        yield i, kmer
