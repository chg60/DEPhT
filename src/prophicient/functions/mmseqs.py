import pathlib

from src.prophicient.functions.run import run


def mmseqs_createdb(fasta, mmseqsdb):
    """
    Creates an MMseqs2 database at `mmseqsdb` using `fasta` as its input.

    :param fasta: path to the FASTA file containing sequences to assemble into an MMseqs2 database
    :type fasta: pathlib.Path
    :param mmseqsdb: path to the desired MMseqs2 database
    :type mmseqsdb: pathlib.Path
    :return:
    """
    c = f"mmseqs createdb {fasta} {mmseqsdb} -v 3"
    run(c)


def mmseqs_cluster(sequence_db, cluster_db, tmp_dir, clustermode, clustersteps,
                   sensitivity, minseqid, coverage, evalue):
    """
    Clusters an MMseqs2 sequence database using the given parameters.

    :param sequence_db: path to the MMseqs2 sequence database to cluster
    :type sequence_db: pathlib.Path
    :param cluster_db: path to the MMseqs2 output cluster database
    :type cluster_db: pathlib.Path
    :param tmp_dir: temporary directory that MMseqs2 can use
    :type tmp_dir: pathlib.Path
    :param clustermode: mmseqs cluster --cluster-mode
    :type clustermode: int
    :param clustersteps: mmseqs cluster --cluster-steps
    :type clustersteps: int
    :param sensitivity: mmseqs cluster -s
    :type sensitivity: float
    :param minseqid: mmseqs cluster --min-seq-id
    :type minseqid: float
    :param coverage: mmseqs cluster -c
    :type coverage: float
    :param evalue: mmseqs cluster -e
    :type evalue: float
    :return:
    """
    c = f"mmseqs cluster {str(sequence_db)} {str(cluster_db)} {str(tmp_dir)} -v 3 " \
        f"--max-seqs 1000 --cluster-mode {clustermode} --cluster-steps {clustersteps} " \
        f"-s {sensitivity} --min-seq-id {minseqid} -c {coverage} -e {evalue}"
    run(c)


def mmseqs_result2profile(sequence_db, cluster_db, profile_db):
    """
    Converts an MMseqs2 cluster output database to a profile database.

    :param sequence_db: path to the MMseqs2 sequence database that was clustered
    :type sequence_db: pathlib.Path
    :param cluster_db: path to the MMseqs2 cluster database
    :type cluster_db: pathlib.Path
    :param profile_db: path to the desired MMseqs2 profile database
    :type profile_db: pathlib.Path
    :return:
    """
    c = f"mmseqs result2profile {str(sequence_db)} {str(sequence_db)} {str(cluster_db)} {str(profile_db)} -v 3"
    run(c)


def mmseqs_profile2consensus(profile_db, consensus_db):
    """
    Extracts consensus sequences from an MMseqs2 profile database and creates
    an MMseqs2 sequence database from those consensus sequences.

    :param profile_db: path to the MMseqs2 profile database
    :type profile_db: pathlib.Path
    :param consensus_db: path to the desired MMseqs2 consensus sequence database
    :type consensus_db: pathlib.Path
    :return:
    """
    c = f"mmseqs profile2consensus {profile_db} {consensus_db} -v 3"
    run(c)


def mmseqs_search(profile_db, consensus_db, align_db, tmp_dir, minseqid, coverage, evalue):
    """
    Searches an MMseqs2 profile database against an MMseqs2 profile consensus
    sequence database, for HMM-based clustering.

    :param profile_db: path to the MMseqs2 profile database
    :type profile_db: pathlib.Path
    :param consensus_db: path to the MMseqs2 consensus sequence database
    :type consensus_db: pathlib.Path
    :param align_db: path to the desired MMseqs2 search result database
    :type align_db: pathlib.Path
    :param tmp_dir: temporary directory that MMseqs2 can use
    :type tmp_dir: pathlib.Path
    :param minseqid: mmseqs search --min-seq-id
    :type minseqid: float
    :param coverage: mmseqs search -c
    :type coverage: float
    :param evalue: mmseqs search --e-profile
    :type evalue: float
    :return:
    """
    c = f"mmseqs search {profile_db} {consensus_db} {align_db} {tmp_dir} -v 3 " \
        f"--max-seqs 1000 --min-seq-id {minseqid} -c {coverage} --cov {coverage} " \
        f"-e {evalue} --e-profile {evalue} --add-self-matches"
    run(c)


def mmseqs_clust(consensus_db, align_db, result_db):
    """
    Clusters the MMseqs2 result database from running 'mmseqs search'.

    :param consensus_db: path to the MMseqs2 consensus sequence database
    :type consensus_db: pathlib.Path
    :param align_db: path to the MMseqs2 search result database
    :type align_db: pathlib.Path
    :param result_db: path to the desired MMseqs2 clustered result database
    :type result_db: pathlib.Path
    :return:
    """
    c = f"mmseqs clust {consensus_db} {align_db} {result_db} -v 3"
    run(c)


def mmseqs_createseqfiledb(sequence_db, cluster_db, sf_db):
    """
    Creates an MMseqs2 seqfileDB from a clustering output database.

    :param sequence_db: the path to the database used as the query in MMseqs2 clustering
    :type sequence_db: pathlib.Path
    :param cluster_db: the path to the MMseqs2 clustering result database
    :type cluster_db: pathlib.Path
    :param sf_db: the path to the desired seqfileDB
    :type sf_db: pathlib.Path
    :return:
    """
    c = f"mmseqs createseqfiledb {sequence_db} {cluster_db} {sf_db} -v 3"
    run(c)


def mmseqs_result2flat(query_db, subject_db, result_db, output):
    """
    Converts an MMseqs2 result seqfileDB to a FASTA-like output file.

    :param query_db: the database used as the query in MMseqs2 clustering
    :type query_db: pathlib.Path
    :param subject_db: the database used as the subject in MMseqs2 clustering
    :type subject_db: pathlib.Path
    :param result_db: the MMseqs2 clustering result database
    :type result_db: pathlib.Path
    :param output: the desired clustering output file
    :type output: pathlib.Path
    :return:
    """
    c = f"mmseqs result2flat {query_db} {subject_db} {result_db} {output} -v 3"
    run(c)


def parse_mmseqs(filepath):
    """
    Parses the indicated MMseqs2 output into a dictionary mapping
    phamids to the geneids that were put in them

    :param filepath: path to the MMseqs2 output file to parse
    :type filepath: pathlib.Path
    :return: phams
    """
    phams = dict()
    phamid = 0
    pham_genes = list()

    with filepath.open("r") as fh:
        prior = fh.readline()
        current = fh.readline()

        # While loop to iterate until EOF
        while current:
            if current.startswith(">"):
                # If current & prior both header lines - new pham block
                if prior.startswith(">"):
                    try:
                        pham_genes.pop(-1)
                    except IndexError:
                        pass
                    phams[phamid] = pham_genes
                    phamid += 1
                    pham_genes = list()
                pham_genes.append(current.lstrip(">").rstrip())
            else:
                # Do nothing on translation lines
                pass
            prior, current = current, fh.readline()
        # Need to dump the last p into the dictionary
        phams[phamid] = pham_genes
    phams.pop(0)  # 0th pham is placeholder

    return phams
