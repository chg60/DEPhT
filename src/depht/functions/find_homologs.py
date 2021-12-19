from depht.classes.hhresult import HHResult
from depht.functions.fasta import write_fasta
from depht.functions.multiprocess import parallelize
from depht.functions.run_command import run_command

HHSEARCH_EVALUE = 1E-04
HHSEARCH_PROB = 90
HHSEARCH_COV = 50


def hhsearch(query, outfile, db, evalue=HHSEARCH_EVALUE, prob=HHSEARCH_PROB,
             cov=HHSEARCH_COV):
    """
    Runs a single instance of hhsearch using the prescribed paths
    and e-value.

    :param query: the filepath to query against an HHSuite3 database
    :type query: pathlib.Path
    :param outfile: the desired output filepath
    :type outfile: pathlib.Path
    :param db: the HHSuite3 database to query
    :type db: pathlib.Path
    :param evalue: the e-value cutoff to use
    :type evalue: float
    :param prob: the probability cutoff to use
    :type prob: float
    :param cov: the coverage cutoff to use
    :type cov: float
    """
    command = f"hhsearch -i {query} -d {db} -o {outfile} -e {evalue} -E " \
              f"{evalue} -p {prob} -cov {cov}"
    run_command(command)


def find_single_homologs(header, sequence, db, tmp_dir, prob=HHSEARCH_PROB,
                         cov=HHSEARCH_COV):
    """
    Runs hhsearch to find functionally annotated homolog(s) for a
    single protein sequence in the indicated database.

    Returns a tuple containing the header and the label of the
    best-scoring match with probability greater than `prob`, if such a
    match exists.

    :param header: a label for the sequence to search
    :type header: str
    :param sequence: the sequence to find homologs for
    :type sequence: str
    :param db: the database to find homologs in
    :type db: pathlib.Path
    :param tmp_dir: the directory where temporary files can go
    :type tmp_dir: pathlib.Path
    :param prob: probability threshold to keep a match
    :type prob: int or float
    :return: hhresult
    """
    # Set up file paths
    query_file = tmp_dir.joinpath(f"{header}.fasta")
    output_file = query_file.with_suffix(".hhr")

    # Write the fasta file, and run hhsearch
    write_fasta([header], [sequence], query_file)
    hhsearch(query_file, output_file, db)

    # Parse the hhsearch result
    hhresult, best_match, best_probability = HHResult(output_file), None, 0.0
    hhresult.parse_result()

    # Cull low-probability matches
    hhresult.matches = [match for match in hhresult.matches if
                        (float(match.probability) >= prob) and
                        ((float(match.match_cols) / float(match.hit_length))
                         * 100 >= cov)]
    if hhresult.matches:
        # Best match is the one with the highest bit-score
        hhresult.matches.sort(key=lambda x: float(x.score), reverse=True)
        best_match = hhresult.matches[0].target_id
        best_probability = hhresult.matches[0].probability
    else:
        if output_file.is_file():
            output_file.unlink()

    query_file.unlink()

    return header, best_match, best_probability


def find_batch_homologs(headers, sequences, db, tmp_dir, cpus):
    """
    Runs hhsearch to find functionally annotated homologs for each
    of the given translations in the indicated database.

    :param headers: labels for the sequence to search
    :type headers: list of str
    :param sequences: the sequences to find homologs for
    :type sequences: list of str
    :param db: the database to find homologs in
    :type db: pathlib.Path
    :param tmp_dir: the directory where temporary files can go
    :type tmp_dir: pathlib.Path
    :param cpus: how many CPU cores to use?
    :type cpus: int
    :return: homologs
    """
    jobs = list()
    for header, sequence in zip(headers, sequences):
        jobs.append((header, sequence, db, tmp_dir))
    homologs = parallelize(jobs, cpus, find_single_homologs)

    return homologs


def find_homologs(contigs, prophage_coords, db, tmp_dir, cpus, min_length=150,
                  cache_scores=True):
    """
    Convenience function for finding all prophage-predicted gene homologs
    across all contigs of a genome.

    Updates contig features in-place.

    :param contigs: the contigs from an input file
    :type contigs: list of Bio.SeqRecord.SeqRecord
    :param prophage_coords: predicted prophage coords on each contig
    :type prophage_coords: list of list of tuple(int, int)
    :param db: the database to find homologs in
    :type db: pathlib.Path
    :param tmp_dir: the directory where temporary files can go
    :type tmp_dir: pathlib.Path
    :param cpus: how many CPU cores to use?
    :type cpus: int
    :param min_length: don't find homologs for genes shorter than this
    :type min_length: int
    :param cache_scores: Toggles whether to store probabilities in the contig
    :type cache_scores: bool
    """
    # Iterate over contigs and prophage coordinate predictions together
    for contig, contig_prophage_coords in zip(contigs, prophage_coords):
        map_geneid_to_feature = dict()
        batch_geneids, batch_sequences = list(), list()
        for i, feature in enumerate(contig.genes):
            geneid = contig.gene_ids[i]

            translation = feature.qualifiers["translation"][0]

            # Only hhsearch if length > min_length
            if len(translation) < min_length:
                continue

            # Only check features that overlap or are in prophages or haven't been
            # searched before
            search = False
            if __feature_in_prophage(feature, contig_prophage_coords):
                search = True
                if contig.hhsearch_scores:
                    if contig.hhsearch_scores[i] > 0:
                        search = False

            if search:
                map_geneid_to_feature[geneid] = feature
                batch_geneids.append(geneid)
                batch_sequences.append(translation)

        homologs = find_batch_homologs(
            batch_geneids, batch_sequences, db, tmp_dir, cpus)

        hhsearch_scores = [0.0] * len(contig.gene_ids)
        for geneid, product, prob in homologs:
            if product:
                feature = map_geneid_to_feature[geneid]
                feature.qualifiers["product"] = [product]

                hhsearch_scores[contig.gene_ids.index(geneid)] = float(prob)

        if cache_scores:
            contig.update_hhsearch_scores(hhsearch_scores)


def __feature_in_prophage(feature, contig_prophage_coords):
    """
    Checks whether the indicated feature overlaps the coordinates of
    a predicted prophage.

    :param feature: the feature to check
    :type feature: Bio.SeqFeature.SeqFeature
    :param contig_prophage_coords: predicted prophage coordinates
    :type contig_prophage_coords: list of tuple(int, int)
    """
    for start, end in contig_prophage_coords:
        if (start < feature.location.start < end) or \
                (start < feature.location.end < end):
            return True
    return False
