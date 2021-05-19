import math

from Bio.SeqFeature import SeqFeature, FeatureLocation
from networkx import DiGraph

from prophicient.classes import kmers
from prophicient.functions.blastn import blastn
from prophicient.functions.fasta import write_fasta


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
KMER_SIZE = 5
FPP = 0.0001

L_SEQ_NAME = "putative_attL_region"
R_SEQ_NAME = "putative_attR_region"

DEFAULTS = {"k": 5, "fpp": 0.0001, "outfmt": 10}


# KMER COUNTING FUNCTIONS
# -----------------------------------------------------------------------------
def find_attachment_site(prophage, l_seq, r_seq, l_origin, r_origin, tmp_dir,
                         k=KMER_SIZE, method="blast", l_name=L_SEQ_NAME,
                         r_name=R_SEQ_NAME):
    """Given the sequences of a putative attL region and putative attR region,
    find the most probable attachment site, dictated by the sequence's length
    and it's distance from the predicted origin position.

    :param prophage: Prophage object to find an attachment site for
    :type prophage: prophicient.classes.prophage.Prophage
    :param l_seq: The sequence of a putative attL region.
    :type l_seq: str
    :param r_seq: The sequence of a putative attR region.
    :type r_seq: str
    :param l_origin: The origin position for the contig's leftmost position
    :type l_origin: int
    :param r_origin: The origin position for the contig's rightmost position
    :type r_origin: int
    :param tmp_dir: The directory for which to write sequences and files to
    :type tmp_dir: pathlib.Path
    :param k: Length of the word size storerd in the DeBruijn graph.
    :type k: int
    :param method: The method to use to find the attachment sequence.
    :type method: str
    :param l_name: Name to give to the putative attL region sequence.
    :type l_name: str
    :param r_name: Name to give to the putative attR region sequence.
    :type r_name: str
    :return: A tuple of information associated with the detected att site.
    :rtype: tuple
    """
    # If method is BLASTn
    if method == "blast":
        # Write the putative attL region to file
        l_seq_path = tmp_dir.joinpath(l_name).with_suffix(".fasta")
        write_fasta([l_name], [l_seq], l_seq_path)

        # Write the putative attR region to file
        r_seq_path = tmp_dir.joinpath(r_name).with_suffix(".fasta")
        write_fasta([r_name], [r_seq], r_seq_path)

        # Use BLASTn to retrieve putative attachment site sequences
        kmer_contigs = blast_attachment_site(
            l_seq_path, r_seq_path, tmp_dir, k=k)
    # If method is DeBruijn graph
    elif method == "graph":
        # Use a DeBruijn to retrieve putative attachment site sequences
        kmer_contigs = graph_attachment_site(l_seq, r_seq, k=k)
    else:
        raise NotImplementedError(
                            f"Attachment site detection method '{method}'"
                            "is not an implemented algorithm")

    if not kmer_contigs:
        return

    # Score putative attachment site sequences
    # TODO paramaterize and allow access to change the 'exponent' variable
    scored_kmer_contigs = [
            (kmer_contig, score_kmer(prophage, kmer_contig, l_origin,
                                     r_origin, k, 3))
            for kmer_contig in kmer_contigs]

    # Sort attachment site sequences by score
    scored_kmer_contigs.sort(key=lambda x: x[1], reverse=True)

    kmer_contig, score = scored_kmer_contigs[0]

    # Create a SeqFeature from the putative attL sequence
    att_l_start = kmer_contig[1] + 1
    att_l_end = kmer_contig[1] + len(kmer_contig[0])
    att_l_feature = SeqFeature(FeatureLocation(att_l_start, att_l_end),
                               strand=1, type="misc_recomb")

    # Create a SeqFeature from the putative attR sequence
    att_r_end = kmer_contig[2]
    att_r_start = kmer_contig[2] - len(kmer_contig[0])
    att_r_feature = SeqFeature(FeatureLocation(att_r_start, att_r_end),
                               strand=1, type="misc_recomb")

    return att_l_feature, att_r_feature, score, kmer_contig[0]


def blast_attachment_site(l_seq_path, r_seq_path, tmp_dir, k=KMER_SIZE):
    """Given the path to files containing the putative attL region and
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
        l_seq_path, r_seq_path, tmp_dir, mode="subject", word_size=k)

    kmer_contigs = []
    for result in blast_results:
        # Append the contig and its positions to the list
        kmer_contig = (result["qseq"],
                       int(result["qstart"]), int(result["send"]))
        kmer_contigs.append(kmer_contig)

    return kmer_contigs


def score_kmer(prophage, kmer_contig, l_origin, r_origin, base, exponent):
    """Score kmer contigs by their relative length and distance from
    an origin position.

    :param kmer_contig: A tuple containing the contig, and its positions
    :type kmer_contig: tuple(str, int, int)
    :param l_origin: The origin position for the contig's leftmost position
    :type l_origin: int
    :param r_origin: The origin position for the contig's rightmost position
    :type r_origin: int
    :param base: The base of the logarithmic function from origin.
    :type base: int
    :param exponent: The exponent of the exponential function from origin
    :type exponent: int
    :return: The score of the given kmer, with respect to the given origins
    :rtype: int
    """
    # Get the distance of the kmer from its left position to the left origin
    l_distance = abs(l_origin - kmer_contig[1])
    # Get the distance of the kmer from its right position to the right origin
    r_distance = abs(r_origin - kmer_contig[2])

    avg_distance = (l_distance + r_distance) / 2

    # Return the length of the kmer, penalized by some function
    # with the average distance as the input
    return len(kmer_contig[0]) - (math.log(avg_distance, base) ** exponent)


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

    # Generator that returns k-mers and their 0-indexed position in the sequence
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        yield i, kmer
