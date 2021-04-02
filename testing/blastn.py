import csv
import shlex
from subprocess import Popen, PIPE


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------

DEFAULT_BLAST_CSV_HEADER = ["qstart", "qend", "sstart", "send",
                            "qseq", "sseq", "length", "mismatch", "gapopen"]

DEFAULT = {"outfmt": 10, "blast_csv_header": DEFAULT_BLAST_CSV_HEADER}


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def locate_subsequence(child_seq_path, parent_seq_path, csv_path):
    """Performs BLASTn on a query and target sequence and identifies
    the coordinates of the section of the parent sequence that corresponds
    to the best alignment.

    :param child_seq_path: Filepath to a fasta-formatted target sequence.
    :type child_seq_path: pathlib.Path
    :param parent_seq_path: Filepath to a fasta-formatted query sequence.
    :type parent_seq_path: pathlib.Path
    :param aln_path: Filepath to dump tabular results of the BLASTn alignment.
    :type aln_path: pathlib.Path
    """
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


def blastn(query, target, out, outfmt=DEFAULT["outfmt"],
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
    command = (f"""blastn -query {query} -subject {target} -out {out} """
               f"""-outfmt "10 {' '.join(header)}" """)

    split_command = shlex.split(command)
    with Popen(args=split_command, stdin=PIPE) as process:
        out, err = process.communicate()


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
