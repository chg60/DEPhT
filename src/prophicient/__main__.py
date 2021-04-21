"""
Prophicient scans bacterial genomes looking for prophages. Regions
identified as prophage candidates are further scrutinized, and
attachment sites identified as accurately as possible before
prophage extraction and generating the final report.
"""
import argparse
import pathlib
import sys
from datetime import date, datetime

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from prophicient.functions import blastn
from prophicient.functions.fasta import write_fasta
from prophicient.functions.att import kmer_count_attachment_site
from prophicient.functions.find_homologs import find_homologs
from prophicient.functions.gene_prediction import *
from prophicient.functions.multiprocess import CPUS, parallelize
# from prophicient.functions.wrapper_basic import autoannotate
from prophicient.functions.prefilter import prefilter_genome, realign_subrecord
from prophicient.functions.prophage_prediction import *
from prophicient.functions.visualization import prophage_diagram

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
TEMP_DIR = pathlib.Path("/tmp/Prophicient/")
PACKAGE_DIR = pathlib.Path(__file__).resolve().parent
DATE = date.today().strftime("%d-%b-%Y").upper()
ANNOTATIONS = {"molecule_type": "DNA", "topology": "linear",
               "data_file_division": "PHG", "date": DATE,
               "accessions": [], "sequence_version": "1",
               "keywords": [], "source": "",
               "organism": "", "taxonomy": [],
               "comment": [""]}
MIN_LENGTH = 20000      # Don't annotate short contigs
META_LENGTH = 100000    # Medium-length contigs -> use metagenomic mode
EXTENSION = 10000


def parse_args(arguments):
    """
    Parse command line arguments.

    :param arguments: command line arguments that program was invoked with
    :type arguments: list
    :return: parsed_args
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", type=pathlib.Path,
                        help="path to a FASTA nucleotide sequence file to "
                             "scan for prophages")
    parser.add_argument("outdir", type=pathlib.Path,
                        help="path where output files should be written")
    parser.add_argument("--gff3", type=pathlib.Path,
                        help="path to a GFF3 file, bypasses auto-annotation")
    parser.add_argument("--verbose", action="store_true",
                        help="toggles verbosity of pipeline")
    parser.add_argument("--no-graphics", action="store_true",
                        help="don't output genome map PDFs for identified "
                             "prophages")
    parser.add_argument("--cpus", type=int, default=CPUS,
                        help=f"number of processors to use [default: {CPUS}]")
    return parser.parse_args(arguments)


def main(arguments):
    """
    Main function that interfaces with command line args and the
    program workflow.

    :param arguments: command line arguments
    :type arguments: list
    """
    args = parse_args(arguments)

    # Verify that the input filepath is valid
    infile = args.infile
    if not infile.is_file():
        print(f"'{str(infile)}' is not a valid input file - exiting...")
        sys.exit(1)

    outdir = args.outdir
    if not outdir.is_dir():
        print(f"'{str(outdir)}' does not exist - creating it...")
        outdir.mkdir(parents=True)

    gff3 = args.gff3
    if gff3 and not gff3.is_file():
        print(f"GFF3 file '{str(gff3)}' does not exist - rollback to "
              f"auto-annotate...")
        gff3 = None

    cpus = args.cpus
    verbose = args.verbose
    diagram = not args.no_graphics

    # Mark program start time
    mark = datetime.now()
    find_prophages(infile, outdir, gff3, cpus, verbose, diagram)
    # execute_prophicient(args.infile, args.outdir, args.database,
    #                     cores=args.cpus, verbose=args.verbose)
    print(f"\nTotal runtime: {str(datetime.now() - mark)}")


def find_prophages(fasta, outdir, gff3=None, cpus=CPUS, verbose=False,
                   diagram=True, extension=EXTENSION):
    """
    Runs through all steps of prophage prediction:

    * auto-annotation with Pyrodigal (skip if gff3)
    * predict prophage genes using binary classifier
    * identify prophage regions
    * identify phage genes in prophage regions
    * detect attL/attR
    * extract final prophage sequences

    :param fasta: the path to a fasta nucleotide sequence file
    containing a mycobacterial genome to find prophages in
    :type fasta: pathlib.Path
    :param outdir: the path to a directory where output files should
    be written
    :type outdir: pathlib.Path
    :param gff3: (optional) the path to a GFF3 file with an annotation
    for the indicated fasta file
    :type gff3: pathlib.Path
    :param cpus: the maximum number of processors to use
    :type cpus: int
    :param verbose: should progress messages be printed along the way?
    :type verbose: bool
    :param diagram: should genome diagrams be created at the end?
    :type diagram: bool
    :param extend_by: number of basepairs to extend predicted prophage regions
    :type extend_by: int
    """
    if verbose:
        print("Loading FASTA file...")
    # Mark FASTA load start
    mark = datetime.now()
    # Parse FASTA file - only keep contigs longer than MIN_LENGTH
    contigs = [x for x in SeqIO.parse(fasta, "fasta") if len(x) >= MIN_LENGTH]
    # Print FASTA load time
    print(f"FASTA load: {str(datetime.now() - mark)}")

    if verbose:
        print("Beginning annotation...")
    # Mark annotation start
    mark = datetime.now()
    # TODO: add mechanism for skipping pyrodigal auto-annotation by
    #  associating records with gff3 file features
    if gff3:
        if verbose:
            print(f"\tUsing gff3 file '{str(gff3)}'...")
        # TODO: check each record's length against MIN_LENGTH - still
        #  skip these to be consistent
        pass
    else:
        if verbose:
            print("\tNo gff3 file - using Pyrodigal and Aragorn...")
        for contig in contigs:
            # Annotate record CDS & t(m)RNA features in-place
            annotate_contig(contig, len(contig) < META_LENGTH)
    # Print annotation time
    print(f"Annotation: {str(datetime.now() - mark)}")

    if verbose:
        print("Looking for high-probability prophage regions...")
    # Mark prophage prediction start
    mark = datetime.now()
    # Get dataframes of CDS features for binary classification
    dataframes = [contig_to_dataframe(contig) for contig in contigs]
    # Perform binary classification of contig CDS features
    gene_predictions = [predict_prophage_genes(df) for df in dataframes]
    # Initial pass at prophage identification
    prophage_predictions = [predict_prophage_coords(x, y, extend_by=extension)
                            for x, y in zip(contigs, gene_predictions)]

    initial_prophages = []
    for contig_index, contig in enumerate(contigs):
        contig_predictions = prophage_predictions[contig_index]
        for prophage_index, prophage_coordinates in enumerate(
                                                    contig_predictions):
            prophage_id = "".join(["prophi", contig.id,
                                   "-", str((prophage_index+1))])
            start = prophage_coordinates[0]
            end = prophage_coordinates[1]

            prophage = Prophage(contig.seq, prophage_id, start=start, end=end)
            initial_prophages.append(prophage)

    # Print prophage prediction time
    print(f"Prediction: {str(datetime.now() - mark)}")

    if len(initial_prophages) == 0:
        print(f"No complete prophages found in {str(fasta)}. PHASTER may "
              f"be able to identify partial (dead) prophages.")
        return

    TEMP_DIR.mkdir(exist_ok=True)

    reference_db_path = pathlib.Path("references/Mycobacteria")
    trimmed_prophages = detect_att_sites(initial_prophages, reference_db_path,
                                         extension, TEMP_DIR)

    final_prophages = list()
    # TODO: add parallel hhsearch to find essential phage functions:
    #  overwrite "hypothetical function" in cds.qualifiers["product"][0]
    # TODO: add attachment core detection
    # TODO: clean up final prophage annotations (add qualifiers for gene and
    #  locus_tag, and maybe gene features for each CDS/tRNA/tmRNA)
    # TODO: store final annotated records in `trimmed_prophages`

    if verbose:
        print("Generating final reports...")
    mark = datetime.now()
    for prophage in final_prophages:
        genbank_filename = outdir.joinpath(f"{prophage.name}.gbk")
        SeqIO.write(prophage, genbank_filename, "genbank")
        if diagram:
            diagram_filename = outdir.joinpath(f"{prophage.name}.pdf")
            prophage_diagram(prophage, diagram_filename)
    print(f"File Dumps: {str(datetime.now() - mark)}")


class Prophage:
    def __init__(self, parent_seq, seq_id, start=None, end=None, strand=1):
        self.parent_seq = parent_seq
        self.id = seq_id

        self.start = None
        self.end = None
        self.strand = strand

        self.feature = None
        self.seq = None
        self.record = None

        self.set_coordinates(start, end)

    def set_coordinates(self, start, end):
        self.start = start
        self.end = end

        self.update_sequence_attributes()

    def set_strand(self, strand):
        self.strand = strand

        self.update_sequence_attributes()

    def update_sequence_attributes(self):
        if self.start is None or self.end is None:
            return

        self.feature = SeqFeature(FeatureLocation(self.start, self.end),
                                  strand=1,
                                  type="source")
        self.seq = self.feature.extract(self.parent_seq)
        self.record = SeqRecord(self.seq, id=self.id)


def execute_prophicient(infile, outdir, database, cores, verbose=False):
    # Parse the input FASTA file
    with infile.open("r") as fasta_reader:
        record_iterator = SeqIO.parse(fasta_reader, "fasta")
        for record in record_iterator:
            record_dir = outdir.joinpath(str(record.name))
            record_dir.mkdir()

            annotate_record(record, record_dir)

            gene_dense_records_and_features = prefilter_genome(record)

            putative_prophage_regions = evaluate_regions_of_interest(
                                            gene_dense_records_and_features,
                                            record_dir, database, cores=cores,
                                            verbose=verbose)

            extracted_prophages = list()
            for region_data, integrase_homolog in putative_prophage_regions:
                prophage_data = extract_prophage(region_data[0],
                                                 integrase_homolog,
                                                 verbose=verbose)

                prophage_global_start = (region_data[1].location.start +
                                         prophage_data[1])
                prophage_global_end = (region_data[1].location.start +
                                       prophage_data[2])

                if verbose:
                    print("Prophage detected at: "
                          f"({prophage_global_start}, {prophage_global_end})")
                    print("...Att sequence detected:\n"
                          f"\t\t{prophage_data[3]}")

                extracted_prophages.append(prophage_data[0])


def extract_prophage(record, integrase_feature, verbose=False):
    l_region = str(record.seq[:integrase_feature.location.start])
    r_region = str(record.seq[integrase_feature.location.start+1:])

    attL, attR = kmer_count_attachment_site(l_region, r_region)
    print(attL)
    att = str(attL.extract(record.seq))

    prophage_start = attL.location.start
    prophage_end = integrase_feature.location.start + attR.location.end
    prophage_feature = SeqFeature(FeatureLocation(prophage_start,
                                                  prophage_end))
    prophage_seq = prophage_feature.extract(record.seq)

    prophage_record = SeqRecord(prophage_seq)
    prophage_record.id = record.id

    realign_subrecord(record, prophage_record, prophage_start, prophage_end)
    prophage_record.features.sort(key=lambda x: x.location.start)

    return prophage_record, prophage_start, prophage_end, att


def detect_att_sites(prophages, reference_db_path, extension,
                     temp_dir):
    for prophage in prophages:
        l_sequence = prophage.seq[:extension]
        r_sequence = prophage

        l_sequence_path = temp_dir.joinpath(f"{prophage.id}_l_region.fasta")
        r_sequence_path = temp_dir.joinpath(f"{prophage.id}_r_region.fasta")

        write_fasta(l_sequence_path, ["l_region"], [str(l_sequence)])
        write_fasta(r_sequence_path, ["r_region"], [str(r_sequence)])

        try:
            l_sequence_blast = blastn.blast_references(
                                l_sequence_path, reference_db_path, temp_dir)
        except blastn.SignificantAlignmentNotFound:
            l_sequence_blast = []

        try:
            r_sequence_blast = blastn.blast_references(
                                r_sequence_path, reference_db_path, temp_dir)
        except blastn.SignificantAlignmentNotFound:
            r_sequence_blast = []
        
        pass
    pass


def evaluate_regions_of_interest(gene_dense_records_and_features, working_dir,
                                 database, cores=1, verbose=False):
    putative_prophage_regions = list()
    for i in range(len(gene_dense_records_and_features)):
        region_tuple = gene_dense_records_and_features[i]

        region_start = region_tuple[1].location.start
        region_end = region_tuple[1].location.end

        if verbose:
            print("Evaluating gene dense region at "
                  f"({region_start}, {region_end})...")

        region_dir = working_dir.joinpath(region_tuple[0].id)
        region_dir.mkdir()

        data_dir = region_dir.joinpath("data")
        data_dir.mkdir()

        integrase_feature = find_integrase_homologs(region_tuple[0], data_dir,
                                                    database, cores=cores,
                                                    verbose=verbose)

        if not integrase_feature:
            continue

        if verbose:
            print("...Integrase homolog found at "
                  f"{(region_start + integrase_feature.location.start)}")

        putative_prophage_regions.append((region_tuple, integrase_feature))

    return putative_prophage_regions


def find_integrase_homologs(record, working_dir, database, cores=1,
                            verbose=False):
    trans_dir = working_dir.joinpath("gene_translations")
    trans_dir.mkdir()

    gene_id_feature_map = dict()
    for feature in record.features:
        if feature.type != "CDS":
            continue

        trans = feature.qualifiers["translation"][0]

        gene_id = feature.qualifiers["locus_tag"][0]
        gene_id_feature_map[gene_id] = feature
        file_path = trans_dir.joinpath(gene_id).with_suffix(".fasta")

        write_fasta(file_path, [gene_id], [trans])

    hhresults_dir = working_dir.joinpath("hhresults")
    hhresults_dir.mkdir()

    integrase_homologs = find_homologs(trans_dir, hhresults_dir, database,
                                       cores=cores, verbose=verbose)

    if not integrase_homologs:
        return

    integrase_homolog = integrase_homologs[0]
    integrase_feature = gene_id_feature_map[integrase_homolog]
    return integrase_feature


def annotate_record(record, working_dir):
    record_file_path = working_dir.joinpath(record.id).with_suffix(".fasta")
    SeqIO.write(record, record_file_path, "fasta")

    data_dir = working_dir.joinpath("data")
    data_dir.mkdir()

    prodigal_out_path = autoannotate(record_file_path, data_dir)
    extract_prodigal_features(record, prodigal_out_path)

    record_file_path = working_dir.joinpath(record.id).with_suffix(".gb")
    record.annotations = ANNOTATIONS
    SeqIO.write(record, record_file_path, "gb")


def extract_prodigal_features(record, prodigal_output_path):
    with prodigal_output_path.open(mode="r") as filehandle:
        prodigal_output = "".join(filehandle.readlines())

    prodigal_genes = PRODIGAL_FORMAT.findall(prodigal_output)

    # Iterate over prodigal_genes (regex "hits" from Prodigal file...)
    for gene in prodigal_genes:
        gene_num = int(gene[0])
        start = int(gene[1]) - 1
        end = int(gene[2])
        strand = int(gene[3])
        partial = int(gene[4])
        start_codon = gene[5]
        rbs_type = gene[6]
        rbs_spacer = gene[7]
        gc_pct = gene[8]

        # Create SeqFeature from these data, and add it to record.features
        qualify = {"note": {"partial": partial, "start_codon": start_codon,
                            "rbs_type": rbs_type, "rbs_spacer": rbs_spacer,
                            "gc content": gc_pct},
                   "product": ["hypothetical protein"],
                   "locus_tag": ["_".join([record.id, "CDS",
                                 str(gene_num+1)])]}
        feature = SeqFeature(FeatureLocation(start, end), strand=strand,
                             qualifiers=qualify, type="CDS")

        feature_sequence = feature.extract(record.seq)
        feature_trans = feature_sequence.translate(to_stop=True, table=11)
        feature.qualifiers["translation"] = [str(feature_trans)]

        record.features.append(feature)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main(sys.argv[1:])
