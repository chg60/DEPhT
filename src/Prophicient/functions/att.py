import csv
import shlex
from subprocess import Popen, PIPE

from Bio.SeqFeature import (FeatureLocation, SeqFeature)
from networkx import DiGraph

from Prophicient.classes import kmers


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT = {"k": 5, "fpp": 0.0001, "outfmt": 10}

BLAST_CSV_HEADER = ["qstart", "qend", "sstart", "send",
                    "qseq", "sseq", "length", "mismatch", "gapopen"]

# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def find_attatchment_site(l_seq, r_seq, working_dir, name=""):
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

    blast_results.sort(key=lambda x: x["length"])
    
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
           header=BLAST_CSV_HEADER):
    command = (f"""blastn -query {query} -subject {target} -out {out} """
               f"""-outfmt "10 {' '.join(BLAST_CSV_HEADER)}" """)
    
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
