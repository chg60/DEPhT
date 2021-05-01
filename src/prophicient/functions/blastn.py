import csv

from prophicient.functions.run import run as run_command


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------

OUTFMT = 10
EVAL_CUTOFF = 0.001
BLAST_HEADER = ["qstart", "qend", "sstart", "send",
                "qseq", "sseq", "sseqid",
                "length", "evalue", "bitscore",
                "gapopen", "mismatch"]

DEFAULTS = {"outfmt": OUTFMT, "blast_csv_header": BLAST_HEADER,
            "eval_cutoff": EVAL_CUTOFF}


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def blast_references(query_seq_path, reference_db_path, temp_dir,
                     eval_cutoff=EVAL_CUTOFF):
    """Performs BLASTn on a query against a reference blast database and
    returns all blast results above a evalue threshold

    :param query_seq_path: Filepath to a fasta-formatted query sequence.
    :type child_seq_path: pathlib.Path
    :param parent_seq_path: Filepath to a BLAST nucleotide database.
    :type parent_seq_path: pathlib.Path
    :param temp_dir: Directory to dump results of the BLASTn alignment.
    :type temp_dir: pathlib.Path
    :param eval_cutoff: Upper E-value threshold to permit results for.
    :type eval_cutoff: float
    :return refined_blast_results: Result dictionaries above E-value cutoff.
    :type refined_blast_results: list[dict]
    """
    csv_path = temp_dir.joinpath(f"{query_seq_path.stem}_blast_results.csv")

    # BLASTn query sequence against reference database
    blastn(query_seq_path, reference_db_path, csv_path, db=True)
    # Read in tabular BLASTn results
    blast_results = read_blast_csv(csv_path)

    if not blast_results:
        raise SignificantAlignmentNotFound(
                    f"BLASTn produced no alignment for {query_seq_path} "
                    f"against {reference_db_path}.")

    refined_blast_results = []
    # Iterate through blast results and remove weaker alignments
    for blast_result in blast_results:
        if float(blast_result["evalue"]) > eval_cutoff:
            continue

        refined_blast_results.append(blast_result)

    return refined_blast_results


def locate_subsequence(child_seq_path, parent_seq_path, temp_dir):
    """Performs BLASTn on a query and target sequence and identifies
    the coordinates of the section of the parent sequence that corresponds
    to the best alignment.

    :param child_seq_path: Filepath to a fasta-formatted target sequence.
    :type child_seq_path: pathlib.Path
    :param parent_seq_path: Filepath to a fasta-formatted query sequence.
    :type parent_seq_path: pathlib.Path
    :param temp_dir: Directory to dump tabular results of the BLASTn alignment.
    :type temp_dir: pathlib.Path
    """
    csv_path = temp_dir.joinpath(f"{child_seq_path.stem}_blast_results.csv")

    # BLASTn child subsequence against parent sequence
    blastn(child_seq_path, parent_seq_path, csv_path)
    # Read in tabular BLASTn results
    blast_results = read_blast_csv(csv_path)

    if not blast_results:
        raise SignificantAlignmentNotFound(
                        f"BLASTn produced no alignment for {parent_seq_path} "
                        f"and {child_seq_path}")

    # Return where the top alignment result aligns to in the parent sequence
    return (int(blast_results[0]["sstart"]), int(blast_results[0]["send"]))


def blastn(query, target, out, db=False, outfmt=OUTFMT,
           header=BLAST_HEADER, word_size=None, gapopen=None, gapextend=None):
    """Performs BLASTn on a query and target sequence.

    :param query: Filepath to a fasta-formatted query sequence.
    :type query: pathlib.Path
    :param target: Filepath to a fasta-formatted target sequence.
    :type target: pathlib.Path
    :param out: Filepath to dump results of the BLASTn alignment.
    :type out: pathlib.Path
    :param outfmt: BLASTn alignment output type
    :type outfmt: int
    :param header: BLASTn tabular results header
    :type header: list[str]
    :param word_size: Word size to use in the BLASTn algorithm
    :type word_size: int
    :param gapopen: Penalty for creating gaps in the alignment
    :type gapopen: int
    :param gapextend: Penalty for extending gaps in the alignment
    :type gapextend: int
    """
    if not db:
        command = (f"""blastn -query {query} -subject {target} -out {out} """
                   f"""-outfmt "10 {' '.join(header)}" """)
    else:
        command = (f"""blastn -query {query} -db {target} -out {out} """
                   f"""-outfmt "10 {' '.join(header)}" """)

    if word_size is not None:
        command = " ".join([command, "-word_size", str(word_size)])
    if gapopen is not None:
        command = " ".join([command, "-gapopen", str(gapopen)])
    if gapextend is not None:
        command = " ".join([command, "-gapextend", str(gapextend)])

    run_command(command)


def read_blast_csv(filepath, header=BLAST_HEADER):
    """Reads in a comma-separated value table.

    :param filepath: Filepath to the tabular results of a BLASTn alignment.
    :type filepath: pathlib.Path
    :param header: BLASTn tabular results header
    :type header: list[str]
    :return: Dictionaries containing key-value pairs, header to value
    :rtype: list[dict]
    """
    blast_results = []
    with filepath.open(mode="r") as filehandle:
        csv_reader = csv.reader(filehandle, delimiter=",", quotechar='"')
        for row in csv_reader:
            row_dict = {header[i]: row[i] for i in range(len(header))}
            blast_results.append(row_dict)

    return blast_results


# ERROR CLASSES
# -----------------------------------------------------------------------------

class SignificantAlignmentNotFound(Exception):
    pass
