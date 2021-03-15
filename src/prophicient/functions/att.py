import csv
import shlex
from subprocess import Popen, PIPE

from Bio.SeqFeature import (FeatureLocation, SeqFeature)
from networkx import DiGraph

from prophicient.classes import kmers


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT = {"k": 5, "fpp": 0.0001, "outfmt": 10}

BLAST_CSV_HEADER = ["qstart", "qend", "sstart", "send",
                    "qseq", "sseq", "length", "mismatch", "gapopen"]


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def blast_attachment_site(l_seq, r_seq, working_dir, name=""):
    l_seq_name = "_".join([name, "attL_region"])
    r_seq_name = "_".join([name, "attR_region"])

    l_seq_path = working_dir.joinpath(l_seq_name).with_suffix(".fasta")
    r_seq_path = working_dir.joinpath(r_seq_name).with_suffix(".fasta")
    outpath = working_dir.joinpath("".join([name, ".csv"]))

    write_fasta(l_seq_path, l_seq, l_seq_name)
    write_fasta(r_seq_path, r_seq, r_seq_name)

    blastn(l_seq_path, r_seq_path, outpath)
    blast_results = read_blast_csv(outpath)

    if not blast_results:
        return None, None

    blast_results.sort(key=lambda x: int(x["length"]), reverse=True)

    att_data = blast_results[0]

    qstart = int(att_data["qstart"]) + 1
    qend = int(att_data["qend"]) + 1

    if qstart > qend:
        temp = qstart
        qstart = qend
        qend = temp
        strand = -1
    else:
        strand = 1

    attL_feature = SeqFeature(FeatureLocation(qstart, qend),
                              strand=strand, type="attL")

    sstart = int(att_data["sstart"]) + 1
    send = int(att_data["send"])

    if sstart > send:
        temp = sstart
        sstart = send
        send = temp
        strand = -1
    else:
        strand = 1

    attR_feature = SeqFeature(FeatureLocation(sstart, send),
                              strand=strand, type="attR")

    if attL_feature.strand != attR_feature.strand:
        return None, None

    return attL_feature, attR_feature


def write_fasta(path, sequence, name):
    with path.open(mode="w") as filehandle:
        filehandle.write("".join([">", name, "\n"]))
        filehandle.write(sequence)


def blastn(query, target, out, outfmt=DEFAULT["outfmt"],
           header=BLAST_CSV_HEADER, k=DEFAULT["k"]):
    command = (f"""blastn -query {query} -subject {target} -out {out} """
               f"""-outfmt "10 {' '.join(BLAST_CSV_HEADER)}" """
               f"-word_size {k}")

    split_command = shlex.split(command)
    with Popen(args=split_command, stdin=PIPE) as process:
        out, err = process.communicate()


def read_blast_csv(filepath, header=BLAST_CSV_HEADER):
    blast_results = []
    with filepath.open(mode="r") as filehandle:
        csv_reader = csv.reader(filehandle, delimiter=",", quotechar='"')
        for row in csv_reader:
            row_dict = {header[i]: row[i] for i in range(len(header))}
            blast_results.append(row_dict)

    return blast_results


def kmer_count_attachment_site(l_seq, r_seq, k=DEFAULT["k"]):
    bfilter = load_bfilter(l_seq, k=k)
    deb_graph = create_debruijn_graph(bfilter, r_seq)
    kmer_contigs = traverse_debruijn_graph(deb_graph, l_seq)

    if not kmer_contigs:
        return None, None

    kmer_contigs.sort(key=lambda x: len(x[0]), reverse=True)

    kmer_contig = kmer_contigs[0]

    attL_start = kmer_contig[1] + 1
    attL_end = kmer_contig[1] + len(kmer_contig)
    attL_feature = SeqFeature(FeatureLocation(attL_start, attL_end),
                              strand=1, type="attL")

    attR_end = kmer_contig[2]
    attR_start = kmer_contig[2] - len(kmer_contig)
    attR_feature = SeqFeature(FeatureLocation(attR_start, attR_end),
                              strand=1, type="attR")

    return attL_feature, attR_feature


def load_bfilter(sequence, k=DEFAULT["k"]):
    bfilter = kmers.BloomFilter(len(sequence))

    for pos, kmer in count_kmers(sequence, k):
        bfilter.add(kmer)

    return bfilter


def create_debruijn_graph(bfilter, sequence, k=DEFAULT["k"]):
    deb_graph = DiGraph()

    anchor_pos = None
    prev_kmer = None
    for pos, kmer in count_kmers(sequence, k):
        if not bfilter.check(kmer):
            anchor_pos = None
            prev_kmer = None
            continue

        if anchor_pos is None:
            anchor_pos = pos

        node = deb_graph.nodes.get(kmer)
        if node is None:
            anchors = set()
            anchors.add(anchor_pos)

            rel_pos_lookup = dict()
            rel_pos = set()
            rel_pos.add(pos - anchor_pos)
            rel_pos_lookup[anchor_pos] = rel_pos

            deb_graph.add_node(kmer, anchors=anchors,
                               rel_pos_lookup=rel_pos_lookup)
        else:
            anchors = node.get("anchors")
            anchors.add(anchor_pos)

            rel_pos_lookup = node.get("rel_pos_lookup")

            rel_pos = rel_pos_lookup.get(anchor_pos, set())
            rel_pos.add(pos - anchor_pos)
            rel_pos_lookup[anchor_pos] = rel_pos

        if prev_kmer is None:
            prev_kmer = kmer
            continue

        edge = deb_graph[prev_kmer].get(kmer)
        if edge is None:
            deb_graph.add_edge(prev_kmer, kmer)

        prev_kmer = kmer

    return deb_graph


def traverse_debruijn_graph(deb_graph, sequence, k=DEFAULT["k"]):
    kmer_contigs = []

    kmer_gen = count_kmers(sequence, k)
    while True:
        try:
            pos, kmer = next(kmer_gen)
        except StopIteration:
            break

        node = deb_graph.nodes.get(kmer)
        if node is not None:
            graph_pos, contig = stitch_kmer_contigs(deb_graph, kmer_gen,
                                                    kmer, k)
            kmer_contigs.append((contig, pos, graph_pos))

    return kmer_contigs


def stitch_kmer_contigs(deb_graph, kmer_gen, kmer, k):
    contig = kmer
    node = deb_graph.nodes[kmer]
    graph_anchor_set = node["anchors"]
    graph_rel_pos_set = build_graph_rel_pos_set(node, graph_anchor_set)
    graph_pos = get_graph_pos(graph_anchor_set, graph_rel_pos_set, k)

    prev_kmer = kmer
    while True:
        try:
            pos, kmer = next(kmer_gen)
        except StopIteration:
            break

        node = deb_graph.nodes.get(kmer)
        if node is None:
            break

        graph_anchor_set = test_valid_kmer_path(deb_graph, graph_anchor_set,
                                                prev_kmer, kmer)
        if not graph_anchor_set:
            break

        graph_rel_pos_set, graph_anchor_set = test_valid_kmer_position(
                                                 deb_graph, graph_rel_pos_set,
                                                 graph_anchor_set, kmer)
        if (not graph_rel_pos_set) or (not graph_anchor_set):
            break

        prev_kmer = kmer
        contig += kmer[-1]
        graph_pos = get_graph_pos(graph_anchor_set, graph_rel_pos_set, k)

    return graph_pos, contig


def test_valid_kmer_path(deb_graph, graph_anchor_set, prev_kmer, kmer):
    node = deb_graph.nodes[kmer]

    edge = deb_graph[prev_kmer].get(kmer)
    if edge is None:
        return

    graph_anchor_set = graph_anchor_set.intersection(node["anchors"])
    return graph_anchor_set


def test_valid_kmer_position(deb_graph, graph_rel_pos_set, graph_anchor_set,
                             kmer):
    inc_graph_rel_pos_set = increment_graph_rel_pos_set(graph_rel_pos_set)

    graph_rel_pos_anchor_set = set()
    int_graph_rel_pos_set = set()

    for anchor, rel_pos_set in deb_graph.nodes[kmer]["rel_pos_lookup"].items():
        for rel_pos in rel_pos_set:
            if rel_pos in inc_graph_rel_pos_set:
                graph_rel_pos_anchor_set.add(anchor)
                int_graph_rel_pos_set.add(rel_pos)

    graph_anchor_set = graph_anchor_set.intersection(graph_rel_pos_anchor_set)

    return int_graph_rel_pos_set, graph_anchor_set


def build_graph_rel_pos_set(node, graph_anchor_set):
    graph_rel_pos_set = set()

    for graph_anchor in graph_anchor_set:
        node_rel_pos_set = node["rel_pos_lookup"][graph_anchor]

        for rel_pos in node_rel_pos_set:
            graph_rel_pos_set.add(rel_pos)

    return graph_rel_pos_set


def increment_graph_rel_pos_set(graph_rel_pos_set):
    inc_graph_rel_pos_set = set()

    for rel_pos in graph_rel_pos_set:
        inc_graph_rel_pos_set.add(rel_pos + 1)

    return inc_graph_rel_pos_set


def get_graph_pos(graph_anchor_set, graph_rel_pos_set, k):
    return max(graph_anchor_set) + max(graph_rel_pos_set) + k


def count_kmers(sequence, k):
    if k >= len(sequence):
        return [(0, sequence)]

    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        yield (i, kmer)
