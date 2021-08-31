import math

from scipy.stats import zscore

from depht.classes.prophage import DEFAULT_PRODUCT
from depht.functions.blastn import blastn, REF_BLASTN_OUTFMT
from depht.functions.fasta import write_fasta
from depht.functions.statistics import transform

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
KMER_SIZE = 5
MIN_ATT_SCORE = 2.2
EVALUE_FILTER = 10000

L_SEQ_NAME = "putative_attL_region"
R_SEQ_NAME = "putative_attR_region"

AQ_WEIGHT = 1
IP_WEIGHT = 0.6
MC_WEIGHT = 0.9
TR_WEIGHT = 0
RC_WEIGHT = 1.5

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

    # Use BLASTn to retrieve putative attachment site sequences
    kmer_contigs = blast_attachment_site(l_seq_path, r_seq_path, tmp_dir, k=k,
                                         evalue=len(l_seq))

    if not kmer_contigs:
        return

    transform_kmer_contig_bitscores(kmer_contigs)

    paired_ref_map = find_reference_att_sites(l_seq_path, r_seq_path,
                                              reference_db_path, tmp_dir,
                                              k, sort_key,
                                              prophage.start, r_seq_start)

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

    attB_table_path = tmp_dir.joinpath("attB.txt")
    dump_reference_attB_sites(paired_ref_map, attB_table_path)

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


def dump_reference_attB_sites(paired_ref_map, outpath):
    filehandle = outpath.open(mode="w")

    for ref_id, ref_data in paired_ref_map.items():
        line_data = [str(ref_data_entry)
                     for ref_data_entry in list(ref_data[4:9]) + [ref_data[3]]]

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
    blast_results = blastn(sequence_path, reference_db_path, tmp_dir,
                           outfmt=REF_BLASTN_OUTFMT)

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
                    att_data = (new_start, new_end, overlap_len, score,
                                ref_id,
                                l_data["sstart"], l_data["send"],
                                r_data["sstart"], r_data["send"])

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
