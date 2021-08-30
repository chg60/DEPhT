import csv

from depht.functions.run_command import run_command

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
BLASTN_OUTFMT = "10 sseqid qstart qend sstart send length gapopen " \
                "mismatch evalue bitscore qseq"

REF_BLASTN_OUTFMT = "10 sseqid qstart qend sstart send length gapopen " \
                    "mismatch evalue bitscore"
BLASTN_EVALUE = 1E-05


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def blastn(query, target, tmp_dir, mode="db", evalue=BLASTN_EVALUE,
           word_size=None, gapopen=None, gapextend=None, outfmt=BLASTN_OUTFMT):
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
    :param word_size: specify a word size (>=4) to use with blastn
    :type word_size: int or None
    :param gapopen: specify a gap open penalty to use with blastn
    :type gapopen: int or None
    :return: results
    """
    # Store any results here
    results = list()

    # Create output filepath
    outfile = tmp_dir.joinpath(f"{query.stem}_blastn_results.csv")

    # Prepare blastn command
    if mode == "db":
        command = f"blastn -query {query} -db {target}"
    elif mode == "subject":
        command = f"blastn -query {query} -subject {target}"
    else:
        raise ValueError("valid blastn modes are: 'db', 'subject'")
    command += f" -out {outfile} -outfmt '{outfmt}'"

    if evalue:
        command += f" -evalue {evalue}"

    if word_size:
        command += f" -word_size {word_size}"

    if gapopen:
        command += f" -gapopen {gapopen}"

    # TODO:
    #  Gap extend and gap open are dependant on each other, for the future,
    #  there needs to be some logic mandating the use of either none or both
    if gapextend:
        command += f" -gapextend {gapextend}"

    run_command(command)

    # Return parsed hits as list of dictionaries
    fields = BLASTN_OUTFMT.split()[1:]
    try:
        blastn_reader = open(outfile, "r")
        csv_reader = csv.DictReader(blastn_reader, fieldnames=fields)
        for row in csv_reader:
            results.append(row)
        blastn_reader.close()
    except FileNotFoundError:
        pass

    return results


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
