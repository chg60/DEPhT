import csv

from prophicient.functions.run_command import run_command


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
BLASTN_OUTFMT = "10 sseqid qstart qend qseq sstart send sseq length gapopen " \
                "mismatch evalue bitscore"
BLASTN_EVALUE = 1E-05


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def blastn(query, target, tmp_dir, mode="db", evalue=BLASTN_EVALUE, **kwargs):
    """
    Runs blastn in either query/subject mode or query/database mode, as
    indicated by `mode`. Returns hits better than `evalue`.

    NOTE: **kwargs will be interpreted as additional blastn parameters.

    :param query: the query FASTA file to use
    :type query: pathlib.Path
    :param target: the target to BLAST against
    :type target: pathlib.Path
    :param tmp_dir: the directory where temporary files can go
    :type tmp_dir: pathlib.Path
    :param mode: how to treat the target (db or subject sequence)
    :type mode: str
    :param evalue: the e-value cutoff to use
    :type evalue: float
    :return:
    """
    # Create output filepath
    outfile = tmp_dir.joinpath(f"{query.stem}_blastn_results.csv")

    # Prepare blastn command
    if mode == "db":
        command = f"blastn -query {query} -db {target}"
    elif mode == "subject":
        command = f"blastn -query {query} -subject {target}"
    else:
        raise ValueError("valid blastn modes are: 'db', 'subject'")
    command += f" -evalue {evalue} -out {outfile} -outfmt '{BLASTN_OUTFMT}'"

    # Interpret any kwargs as blastn keywords
    for key, value in kwargs.items():
        command += f" -{key} {value}"
    run_command(command)

    # Return parsed hits as list of dictionaries
    fields = BLASTN_OUTFMT.split()[1:]
    try:
        blastn_reader = csv.DictReader(open(outfile, "r"), fieldnames=fields)
        return [row for row in blastn_reader]
    except FileNotFoundError:
        return []


def locate_subsequence(query, subject, tmp_dir):
    """
    Uses `blastn` method to find the coordinates of a query sequence
     in a target sequence.

    :param query: the query sequence to search
    :type query: pathlib.Path
    :param subject: the subject sequence to search
    :type subject: pathlib.Path
    :param tmp_dir: the directory where temporary files can go
    :type tmp_dir: pathlib.Path
    :return: start, end
    """
    blast_results = blastn(query, subject, tmp_dir, mode="subject")

    # Return best alignment result location if any alignments were found
    if blast_results:
        start = int(blast_results[0]["sstart"])
        end = int(blast_results[0]["send"])
        return start, end
