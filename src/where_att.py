from pathlib import Path
import time

from Bio import SeqIO
from Bio.SeqFeature import (FeatureLocation, SeqFeature)
from networkx import DiGraph

from prophicient.classes import kmers


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT = {"k": 5, "fpp": 0.0001}

INPUT_PATHS = ["prophiGD17-1.fasta", "prophiGD17-2.fasta",
               "prophiGD12-2.fasta"]

# MAIN FUNCTIONS
# -----------------------------------------------------------------------------


def main(k=DEFAULT["k"]):
    for input_path in INPUT_PATHS:
        input_path = Path.cwd().joinpath(input_path)
        print(f"Finding attatchment site for {input_path.stem}...")
        record = SeqIO.read(input_path, "fasta")

        l_region = str(record.seq[:20000])
        r_region = str(record.seq[-20000:])

        start = time.time()
        attP = find_naive_attatchment_site(l_region, r_region, k=k)
        stop = time.time()

        print("Time elapsed: {:.3f}s".format(stop - start))
        print(f"Longest naive contig @ pos {attP.location.start}: "
              f"{attP.extract(l_region)}")

        start = time.time()
        attP = find_attatchment_site(l_region, r_region, k=k)
        stop = time.time()

        if attP is None:
            print("Traversing algorithm didn't find any contigs")
            continue

        print("Time elapsed: {:.3f}s".format(stop - start))
        print(f"Longest traced contig @ pos {attP.location.start}: "
              f"{attP.extract(l_region)}")

        print("")


def find_attatchment_site(l_sequence, r_sequence, k=DEFAULT["k"]):
    bfilter = load_bloom_filter(l_sequence, k=k)
    deb_graph = create_debruijn_graph(bfilter, r_sequence, k=k)
    contigs = traverse_debruijn_graph(deb_graph, l_sequence, k=k)

    contigs.sort(key=lambda x: len(x[0]), reverse=True)

    if not contigs:
        return None

    att_end = (contigs[0][1] + k)
    att_start = (att_end - len(contigs[0][0]))
    att_feature = SeqFeature(FeatureLocation(att_start, att_end), strand=1,
                             type="attP")

    return att_feature


def find_naive_attatchment_site(l_sequence, r_sequence, k=DEFAULT["k"]):
    bfilter = load_bloom_filter(l_sequence, k=k)
    deb_graph = create_debruijn_graph(bfilter, r_sequence, k=k)
    contigs = trace_debruijn_graph(deb_graph, l_sequence, k=k)

    contigs.sort(key=lambda x: x[1])
    contigs.sort(key=lambda x: len(x[0]), reverse=True)

    att_end = (contigs[0][1] + k)
    att_start = (att_end - len(contigs[0][0]))
    att_feature = SeqFeature(FeatureLocation(att_start, att_end), strand=1,
                             type="attP")

    return att_feature


def trace_debruijn_graph(deb_graph, sequence, k=DEFAULT["k"]):
    contigs = []

    contig = None
    for kmer, pos in count_kmers(sequence, k):
        node = deb_graph.nodes.get(kmer)

        if node is None:
            if contig is not None:
                contigs.append((contig, pos-1))

            contig = None
            continue

        if contig is not None:
            contig = contig + kmer[-1]
        else:
            contig = kmer

    return contigs


def create_debruijn_graph(bfilter, sequence, k=DEFAULT["k"]):
    deb_graph = DiGraph()

    prev_match = False
    prev_kmer = None
    for kmer, pos in count_kmers(sequence, k):
        match = bfilter.check(kmer)

        if match:
            node = deb_graph.nodes.get(kmer)

            if node is None:
                positions = set()
                positions.add(pos)

                deb_graph.add_node(kmer, positions=positions)
            else:
                node["positions"].add(pos)

            if prev_match:
                edge = deb_graph[prev_kmer].get(kmer)

                if edge is None:
                    deb_graph.add_edge(prev_kmer, kmer,
                                       pos_pairs=[(pos-1, pos)])
                else:
                    pos_pairs = edge["pos_pairs"]
                    pos_pairs.append((pos-1, pos))
                    edge["pos_pairs"] = pos_pairs

        prev_match = match
        prev_kmer = kmer

    return deb_graph


def load_bloom_filter(sequence, k=DEFAULT["k"], fpp=DEFAULT["fpp"]):
    bfilter = kmers.BloomFilter(len(sequence), fpp=fpp)

    for kmer, pos in count_kmers(sequence, k):
        bfilter.add(kmer)

    return bfilter


def count_kmers(sequence, k):
    num_kmers = len(sequence) - k
    if num_kmers <= 0:
        raise

    for i in range(num_kmers):
        yield sequence[i:i+k], i


def traverse_debruijn_graph(deb_graph, sequence, k=DEFAULT["k"]):
    contigs = []

    kmers = [kmer[0] for kmer in count_kmers(sequence, k)]
    num_kmers = len(kmers)

    contig = None
    counter = 0
    while True:
        kmer = kmers[counter]
        node = deb_graph.nodes.get(kmer)

        if node is not None:
            contig_kmers, kmer_start = stitch_kmer_path(
                                                deb_graph, kmers[counter:])

            contig = contig_kmers[0]
            if len(contig_kmers) > 1:
                for contig_kmer in contig_kmers[1:]:
                    contig += contig_kmer[-1]

            counter += len(contig_kmers)
            contigs.append((contig, counter-1))
        else:
            counter += 1

        if counter >= num_kmers:
            break

    return contigs


def stitch_kmer_path(deb_graph, kmers):
    contig_kmers = list()

    prev_kmer_start = None
    for i in range(len(kmers)):
        kmer = kmers[i]

        node = deb_graph.nodes.get(kmer)

        if node is None:
            break

        kmer_start = validate_kmer_path(deb_graph, (contig_kmers+[kmer]))
        if kmer_start is None:
            break
        else:
            prev_kmer_start = kmer_start

        contig_kmers.append(kmer)

    return contig_kmers, prev_kmer_start


def validate_kmer_path(deb_graph, kmer_path):
    pos_lookup = create_kmer_position_lookup(deb_graph, kmer_path)

    if not pos_lookup:
        return None

    anchor_kmer = kmer_path[0]
    anchor_node = deb_graph.nodes.get(anchor_kmer)

    valid = True
    for anchor_pos in anchor_node["positions"]:
        valid = True
        for i in range(len(kmer_path)):
            curr_kmer = pos_lookup.get(anchor_pos+i)

            if curr_kmer is None:
                valid = False
                break

            if curr_kmer != kmer_path[i]:
                valid = False
                break

        if valid:
            break

    if valid:
        return anchor_pos


def create_kmer_position_lookup(deb_graph, kmers):
    pos_lookup = dict()

    for kmer in kmers:
        node = deb_graph.nodes.get(kmer)

        if node is None:
            raise

        for pos in node["positions"]:
            pos_lookup[pos] = kmer

    return pos_lookup


if __name__ == "__main__":
    main(k=7)
