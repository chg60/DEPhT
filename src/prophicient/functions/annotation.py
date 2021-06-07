"""
This file contains functions for predicting protein coding and
tRNA/tmRNA genes on bacterial contigs.
"""

import pathlib
from tempfile import mkstemp

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from prophicient.functions.fasta import parse_fasta
from prophicient.functions.run_command import run_command

MIN_LENGTH = 20000      # Don't annotate short contigs
MIN_CDS_FEATURES = 51
META_LENGTH = 100000    # Medium-length contigs -> use metagenomic mode

DEFAULT_PRODUCT = "hypothetical_protein"


def prodigal(infile, outfile, meta=False):
    """
    Runs Prodigal (v2.6.3) on the indicated input file, in metagenomic
    mode if `meta` == True, with predicted CDS features written to the
    indicated output file.

    :param infile: the input file with nucleotide sequence to annotate
    :type infile: pathlib.Path
    :param outfile: the output file where predicted genes should go
    :type outfile: pathlib.Path
    :param meta: run in metagenomic mode?
    :type meta: bool
    """
    command = f"prodigal -a {outfile} -i {infile}"
    if meta:
        command += " -p meta"
    run_command(command)


def aragorn(infile, outfile):
    """
    Runs Aragorn (v1.2.38) on the indicated input file, with predicted
    tRNA and tmRNA features written to the indicated output file.

    :param infile: the input file with nucleotide sequence to annotate
    :type infile: pathlib.Path
    :param outfile: the output file where predicted t(m)RNAs should go
    :type outfile: pathlib.Path
    """
    command = f"aragorn -gcbact -c -d -wa -o {outfile} {infile}"
    run_command(command)


def prodigal_reader(filepath):
    """
    Generator that creates SeqFeatures from genes in the indicated
    Prodigal translation FASTA file.

    :param filepath: the path to a Prodigal translation FASTA file
    :type filepath: pathlib.Path
    """
    headers, sequences = parse_fasta(filepath)
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

        yield ftr


def aragorn_reader(filepath):
    """
    Generator that creates SeqFeatures from rows in the indicated
    Aragorn output file.

    :param filepath: the path to an Aragorn output file
    :type filepath: pathlib.Path
    """
    with open(filepath, "r") as ar:
        next(ar), next(ar)
        for row in ar:
            row = row.rstrip().split()[1:]  # tokenize line, skip index

            # Get coordinates
            if row[1].startswith("c"):
                strand = -1     # reverse oriented
                coords = row[1][2:-1].split(",")
            else:
                strand = 1      # forward oriented
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

            yield ftr


def annotate_contig(contig, tmp_dir, no_trna=False):
    """
    Uses Prodigal to predict protein-coding genes, and Aragorn to
    predict t(m)RNA genes on bacterial contigs. All resultant features
    are appended directly to the contig's features list.

    :param contig: the nucleotide sequence to predict genes on
    :type contig: Bio.SeqRecord.SeqRecord
    :param tmp_dir: temporary directory where files can go
    :type tmp_dir: pathlib.Path
    :param no_trna: don't annotate tRNAs
    :type no_trna: bool
    """
    # Set up to run Prodigal
    infile = mkstemp(suffix=".fna", prefix=f"{contig.id}_", dir=tmp_dir)[-1]
    infile = pathlib.Path(infile)

    with infile.open("w") as prodigal_writer:
        SeqIO.write(contig, prodigal_writer, "fasta")

    prodigal_out = infile.with_suffix(".faa")

    # Run Prodigal
    meta = len(contig) < META_LENGTH
    prodigal(infile, prodigal_out, meta)

    for ftr in prodigal_reader(prodigal_out):
        contig.features.append(ftr)

    if not no_trna:
        # Set up to run Aragorn
        aragorn_out = infile.with_suffix(".txt")

        # Run Aragorn
        aragorn(infile, aragorn_out)

        for ftr in aragorn_reader(aragorn_out):
            contig.features.append(ftr)
