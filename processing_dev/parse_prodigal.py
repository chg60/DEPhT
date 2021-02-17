from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re


# Compile the prodigal pattern
PATTERN = re.compile(">\w+_(\d+) # (\d+) # (\d+) # (-?1) # ID=\w+;partial=(\d+);start_type=(\w+);rbs_motif=(.*|None);rbs_spacer=(.*bp|None);gc_cont=(\d.\d+)\s+")


# Store the genome nucleotide sequence as a SeqRecord
with open("Phoebe.fasta", "r") as fh:
    record = SeqIO.read(fh, "fasta")

# Read the Prodigal output and parse it using PATTERN
with open("Phoebe.gbk", "r") as fh:
    contents = "".join(fh.readlines())
prodigal_genes = PATTERN.findall(contents)

# Iterate over prodigal_genes (regex "hits" from Prodigal file...)
for gene in prodigal_genes:
    gene_num = int(gene[0])
    start = int(gene[1])
    end = int(gene[2])
    strand = int(gene[3])
    partial = int(gene[4])
    start_codon = gene[5]
    rbs_type = gene[6]
    rbs_spacer = gene[7]
    gc_pct = gene[8]

    # Create SeqFeature from these data, and add it to record.features
    qualify={"note":{"partial": partial, "start_codon": start_codon,
            "rbs_type": rbs_type, "rbs_spacer": rbs_spacer,
            "gc content": gc_pct} }
    feature = SeqFeature(FeatureLocation(start, end),strand, qualifiers=qualify)
    record.features.append(feature)
