"""
This file contains functions for predicting protein coding and
tRNA/tmRNA genes on bacterial contigs.
"""

import pathlib
import shlex
import time
from subprocess import Popen, DEVNULL

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from depht.functions.fasta import parse_fasta

MIN_LENGTH = 20000      # Don't annotate short contigs
META_LENGTH = 100000    # Medium-length contigs -> use metagenomic mode

DEFAULT_PRODUCT = "hypothetical protein"


def prodigal(infile, outfile, meta=False):
    """
    Runs Prodigal (v2.6.3) asynchronously on the indicated input file,
    in metagenomic mode if `meta` == True, with predicted CDS features
    written to the indicated output file.

    :param infile: the input file with nucleotide sequence to annotate
    :type infile: pathlib.Path
    :param outfile: the output file where predicted genes should go
    :type outfile: pathlib.Path
    :param meta: run in metagenomic mode?
    :type meta: bool
    :return: process_handle
    """
    try:
        command = f"prodigal -a {outfile} -i {infile} -n -c"
        if meta:
            command += " -p meta"
        command = shlex.split(command)
        process_handle = Popen(args=command, stdout=DEVNULL, stderr=DEVNULL,
                               close_fds=True)
    except OSError:
        raise RuntimeError("Unable to locate Prodigal")
    return process_handle


def parse_prodigal(outfile):
    """
    Parses Prodigal output file into a list of BioPython SeqFeatures.

    :param outfile: the path to the output file written by Prodigal
    :type outfile: pathlib.Path
    :return: features
    """
    features = list()

    headers, sequences = parse_fasta(outfile)
    for header, sequence in zip(headers, sequences):
        header = header.split(" # ")
        start, end, strand = int(header[1]), int(header[2]), int(header[3])
        notes = header[-1].split(";")
        motif = notes[-3].split("=")[-1]
        spacer = notes[-2].split("=")[-1]

        ftr = SeqFeature(location=FeatureLocation(start - 1, end),
                         type="CDS", strand=strand)
        ftr.qualifiers["gene"] = [""]
        ftr.qualifiers["locus_tag"] = [""]
        ftr.qualifiers["note"] = [f"rbs_motif: {motif}; rbs_spacer: {spacer}"]
        ftr.qualifiers["transl_table"] = [11]
        ftr.qualifiers["product"] = [DEFAULT_PRODUCT]
        ftr.qualifiers["translation"] = [sequence.rstrip("*")]

        features.append(ftr)

    return features


def aragorn(infile, outfile):
    """
    Runs Aragorn (v1.2.38) asynchronously on the indicated input file,
    with predicted tRNA and tmRNA features written to the indicated
    output file.

    :param infile: the input file to be used by Aragorn
    :type infile: pathlib.Path
    :param outfile: the output file to be written by Aragorn
    :type outfile: pathlib.Path
    :return: process_handle
    """
    try:
        command = f"aragorn -gcbact -l -d -wa -o {outfile} {infile}"
        command = shlex.split(command)
        process_handle = Popen(args=command, stdout=DEVNULL, stderr=DEVNULL,
                               close_fds=True)
    except OSError:
        raise RuntimeError("Unable to locate Aragorn")
    return process_handle


def parse_aragorn(outfile):
    """
    Parses Aragorn output file into a list of BioPython SeqFeatures.

    :param outfile: the path to the output file written by Aragorn
    :type outfile: pathlib.Path
    :return: features
    """
    features = list()

    aragorn_reader = open(outfile, "r")

    # Skip the 2 header lines
    for _ in range(2):
        next(aragorn_reader)

    # Now parse the remaining lines into Bio.SeqFeature.SeqFeature objects
    for row in aragorn_reader:
        # Tokenize the line; skip the tRNA number, which is meaningless
        row = row.rstrip().split()[1:]

        # Get coordinates
        if row[1].startswith("c"):
            strand = -1  # reverse oriented
            coords = row[1][2:-1].split(",")
        else:
            strand = 1  # forward oriented
            coords = row[1][1:-1].split(",")
        start, end = int(coords[0]), int(coords[1])

        # Check if this is a tRNA or tmRNA
        if row[0] == "tmRNA":
            ftr = SeqFeature(location=FeatureLocation(start - 1, end),
                             type="tmRNA", strand=strand)
            ftr.qualifiers["gene"] = [""]
            ftr.qualifiers["locus_tag"] = [""]
            tag_peptide = row[-1].rstrip("*")
            ftr.qualifiers["note"] = [f"tag peptide: {tag_peptide}"]
        else:
            ftr = SeqFeature(location=FeatureLocation(start - 1, end),
                             type="tRNA", strand=strand)
            ftr.qualifiers["gene"] = [""]
            ftr.qualifiers["locus_tag"] = [""]
            ftr.qualifiers["note"] = [f"{row[0]}{row[-1]}"]
            if "?" in row[0] or "SeC" in row[0] or "Pyl" in row[0]:
                ftr.qualifiers["product"] = ["tRNA-OTHER"]
            else:
                ftr.qualifiers["product"] = [f"{row[0]}"]

        features.append(ftr)

    # Close the file handle
    aragorn_reader.close()

    return features


def annotate_record(record, tmp_dir, trna=True):
    """
    Uses Prodigal to predict protein-coding genes, and Aragorn to
    predict t(m)RNA genes on bacterial contigs. All resultant features
    are appended directly to the contig's features list.

    :param record: the nucleotide sequence to predict genes on
    :type record: Bio.SeqRecord.SeqRecord
    :param tmp_dir: temporary directory where files can go
    :type tmp_dir: pathlib.Path
    :param trna: don't annotate tRNAs
    :type trna: bool
    """
    # Set up to run Prodigal and Aragorn
    infile = tmp_dir.joinpath(f"{record.id}.fna")

    infile_writer = infile.open("w")
    SeqIO.write(record, infile_writer, "fasta")
    infile_writer.close()

    # Name the output files
    aragorn_out = infile.with_suffix(".txt")
    prodigal_out = infile.with_suffix(".faa")

    # Set up to run Prodigal first, since it takes the longest to run
    prodigal_process = prodigal(infile, prodigal_out,
                                len(record) < META_LENGTH)

    # Now kick off Aragorn - wait until it finishes, then begin parsing output
    if trna:
        aragorn_process = aragorn(infile, aragorn_out)
        while aragorn_process.poll() is None:
            time.sleep(0.5)
        for ftr in parse_aragorn(aragorn_out):
            record.features.append(ftr)

    # Now wait until Prodigal finishes to parse its output
    while prodigal_process.poll() is None:
        time.sleep(0.5)
    for ftr in parse_prodigal(prodigal_out):
        record.features.append(ftr)

    # Sort contig features on start position
    record.features.sort(key=lambda x: x.location.start)
