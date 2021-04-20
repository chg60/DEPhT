"""
This file contains functions for predicting protein coding and
tRNA/tmRNA genes on bacterial contigs.
"""

import pathlib

import pyrodigal
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from prophicient.functions.run import run


def annotate_contig(contig, meta=False):
    """
    Annotates CDS and t(m)RNA features on the given contig, in-place.
    Identified features are sorted on their start coordinates.

    :param contig: the nucleotide sequence to predict genes on
    :type contig: Bio.SeqRecord.SeqRecord
    :param meta: run Pyrodigal in metagenomic mode?
    :type meta: bool
    """
    predict_cds_features(contig, meta)
    predict_trna_features(contig)
    contig.features.sort(key=lambda x: x.location.start)


def predict_cds_features(contig, meta=False):
    """
    Uses Pyrodigal (Cythonized Prodigal) to predict protein-coding
    genes on bacterial contigs. Resultant CDS features are appended
    directly to the contig's feature list.

    NOTE: meta=True is slow on full-length bacterial genomes, so only
    use for short contigs (20-100 kb).

    :param contig: the nucleotide sequence to predict genes on
    :type contig: Bio.SeqRecord.SeqRecord
    :param meta: run Pyrodigal in metagenomic mode?
    :type meta: bool
    """
    # Pyrodigal requires string as input
    sequence = str(contig.seq)

    # Create Pyrodigal gene caller, set metagenomic flag as appropriate
    cds_caller = pyrodigal.Pyrodigal(meta=meta)
    if not meta:
        cds_caller.train(sequence)

    for ftr in pyrodigal_output_parser(cds_caller.find_genes(sequence)):
        contig.features.append(ftr)


def predict_trna_features(contig):
    """
    Uses Aragorn (version 1.2.38) to predict tRNA and tmRNA genes on
    bacterial contigs. Resultant t(m)RNA genes are appended directly
    to the contig's feature list.

    :param contig: the nucleotide sequence to predict genes on
    :type contig: Bio.SeqRecord.SeqRecord
    """
    # Set up to run Aragorn
    infile = pathlib.Path("/tmp/aragorn_input.fasta")
    with infile.open("w") as aragorn_writer:
        SeqIO.write(contig, aragorn_writer, "fasta")

    outfile = pathlib.Path("/tmp/aragorn_output.txt")

    # Run Aragorn
    command = f"aragorn -gcbact -c -d -wa -o {str(outfile)} {str(infile)}"
    run(command, verbose=False)

    # Parse Aragorn tRNAs into new SeqFeatures
    for ftr in aragorn_output_parser(outfile):
        contig.features.append(ftr)


def pyrodigal_output_parser(genes):
    """
    Generator that creates SeqFeatures from genes in the given
    Pyrodigal gene predictions.

    :param genes: the genes predicted by pyrodigal
    :type genes: pyrodigal.Genes
    """
    for gene in genes:
        ftr = SeqFeature(location=FeatureLocation(gene.begin - 1, gene.end),
                         type="CDS", strand=gene.strand)
        ftr.qualifiers["product"] = ["hypothetical protein"]
        ftr.qualifiers["note"] = [f"rbs_motif: {gene.rbs_motif}; "
                                  f"rbs_spacer: {gene.rbs_spacer}"]
        ftr.qualifiers["translation"] = [gene.translate(11).rstrip("*")]
        yield ftr


def aragorn_output_parser(filepath):
    """
    Generator that creates SeqFeatures from rows in the indicated
    Aragorn output file.

    :param filepath: the path to an Aragorn output file
    :type filepath: pathlib.Path
    """
    with filepath.open("r") as aragorn_reader:
        next(aragorn_reader), next(aragorn_reader)
        for row in aragorn_reader:
            row = row.rstrip().split()[1:]  # tokenize line, skip index

            # Get coordinates
            if row[1].startswith("c"):
                strand = -1  # Reverse oriented
                coords = row[1][2:-1].split(",")
            else:
                strand = 1  # Forward oriented
                coords = row[1][1:-1].split(",")
            start, end = int(coords[0]), int(coords[1])

            # Check if it's a tRNA or tmRNA
            if row[0] == "tmRNA":
                ftr = SeqFeature(location=FeatureLocation(start - 1, end),
                                 type="tmRNA", strand=strand)
                tag_peptide = row[-1].rstrip("*")
                ftr.qualifiers["note"] = [f"tag peptide: {tag_peptide}"]
            else:
                ftr = SeqFeature(location=FeatureLocation(start - 1, end),
                                 type="tRNA", strand=strand)
                if "?" in row[0] or "SeC" in row[0] or "Pyl" in row[0]:
                    ftr.qualifiers["product"] = ["tRNA-OTHER"]
                else:
                    ftr.qualifiers["product"] = [f"{row[0]}"]
                ftr.qualifiers["note"] = [f"{row[0]}{row[-1]}"]

            yield ftr
