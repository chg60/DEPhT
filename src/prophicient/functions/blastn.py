import csv

from prophicient.functions.run import run as run_command


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------

DEFAULT_BLAST_CSV_HEADER = ["qstart", "qend", "sstart", "send",
                            "qseq", "sseq", "sseqid",
                            "length", "evalue", "bitscore"]

DEFAULT = {"outfmt": 10, "blast_csv_header": DEFAULT_BLAST_CSV_HEADER,
           "eval_cutoff": 0.001}


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def blast_references(query_seq_path, reference_db_path, temp_dir,
                     eval_cutoff=DEFAULT["eval_cutoff"]):
    """Performs BLASTn on a query against a reference blast database and
    returns all blast results above a evalue threshold
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


def blastn(query, target, out, db=False, outfmt=DEFAULT["outfmt"],
           header=DEFAULT["blast_csv_header"]):
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
    """
    if not db:
        command = (f"""blastn -query {query} -subject {target} -out {out} """
                   f"""-outfmt "10 {' '.join(header)}" """)
    else:
        command = (f"""blastn -query {query} -db {target} -out {out} """
                   f"""-outfmt "10 {' '.join(header)}" """)

    run_command(command)


def read_blast_csv(filepath, header=DEFAULT["blast_csv_header"]):
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
