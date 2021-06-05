"""
Prophicient scans bacterial genomes looking for prophages. Regions
identified as prophage candidates are further scrutinized, and
attachment sites identified as accurately as possible before
prophage extraction and generating the final report.
"""
import argparse
import pathlib
import shutil
import sys
from datetime import datetime

from Bio import SeqIO

from prophicient import PACKAGE_DIR
from prophicient.classes.prophage import ANNOTATIONS, Prophage
from prophicient.functions.annotation import annotate_contig, MIN_LENGTH
from prophicient.functions.att import find_attachment_site
from prophicient.functions.blastn import blastn
from prophicient.functions.fasta import write_fasta
from prophicient.functions.find_homologs import find_homologs
from prophicient.functions.mmseqs import assemble_bacterial_mask
from prophicient.functions.multiprocess import PHYSICAL_CORES
from prophicient.functions.prophage_prediction import predict_prophage_coords
from prophicient.functions.visualization import prophage_diagram

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
RUN_MODE_MAP = {"fast": 0, "normal": 1, "strict": 2}
DEFAULT_RUN_MODE = "normal"

TMP_DIR = pathlib.Path("/tmp/prophicient")

PROPHAGE_PREFIX = "prophi"
PROPHAGE_DELIMITER = "-"

EXTEND_BY = 4000
SEARCH_EXTEND_MULT = 2
MIN_SIZE = 10000
REF_BLAST_SORT_KEY = "bitscore"

NORMAL_PRODUCT_THRESHOLD = 3
STRICT_PRODUCT_THRESHOLD = 10


# DEFAULT FILE PATHS
# -----------------------------------------------------------------------------
BLASTN_DB = PACKAGE_DIR.joinpath("data/blastn/mycobacteria")

ESSENTIAL_DB = PACKAGE_DIR.joinpath("data/hhsearch/essential/essential")
EXTENDED_DB = PACKAGE_DIR.joinpath("data/hhsearch/extended/extended")

BACTERIAL_REF_FASTA = PACKAGE_DIR.joinpath("data/mmseqs/bacterial_genes.fasta")
BACTERIAL_REF_VALUES = PACKAGE_DIR.joinpath("data/mmseqs/bacterial_genes.pbv")


def parse_args(arguments):
    """
    Parse command line arguments.

    :param arguments: command line arguments that program was invoked with
    :type arguments: list
    :return: parsed_args
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", type=pathlib.Path,
                        help="Path to a FASTA nucleotide sequence file to "
                             "scan for prophages")
    parser.add_argument("outdir", type=pathlib.Path,
                        help="Path where output files should be written")

    parser.add_argument("--verbose", action="store_true",
                        help="Toggle on verbosity of pipeline")
    parser.add_argument("--no-draw", action="store_false", dest="draw",
                        help="Toggle off genome map PDFs for identified "
                             "prophages")

    parser.add_argument("--mode", type=str,
                        choices=RUN_MODE_MAP.keys(), default=DEFAULT_RUN_MODE,
                        help="Choose run mode of the pipeline.")
    parser.add_argument("--cpus", type=int, default=PHYSICAL_CORES,
                        help=f"Number of processors to use [default: "
                             f"{PHYSICAL_CORES}]")

    parser.add_argument("--dump_data", action="store_true",
                        help="Choose whether to make data accessible")

    parser.add_argument("--prophage-prefix", type=str,
                        default=PROPHAGE_PREFIX,
                        help="Prefix to append to the beginning of the names "
                             "given to detected prophages")
    parser.add_argument("--prophage-delimiter", type=str,
                        default=PROPHAGE_DELIMITER,
                        help="Suffix delimiter to use to separate the names "
                             "given to detected prophages and their number")
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
        return

    outdir = args.outdir
    if not outdir.is_dir():
        print(f"'{str(outdir)}' does not exist - creating it...")
        outdir.mkdir(parents=True)

    if not args.dump_data:
        temp_dir = TMP_DIR
    else:
        temp_dir = outdir.joinpath("data")

    # Create temporary dir, if it doesn't exist.
    if not temp_dir.is_dir():
        temp_dir.mkdir()

    # Mark program start time
    mark = datetime.now()
    find_prophages(infile, outdir, temp_dir,
                   cpus=args.cpus, verbose=args.verbose, draw=args.draw,
                   extend_by=EXTEND_BY, run_mode=RUN_MODE_MAP[args.mode],
                   prefix=args.prophage_prefix,
                   delimiter=args.prophage_delimiter)
    print(f"\nTotal runtime: {str(datetime.now() - mark)}")

    if not args.dump_data:
        # Clean up after ourselves
        shutil.rmtree(TMP_DIR)


def find_prophages(fasta, outdir, tmp_dir, cpus=PHYSICAL_CORES,
                   verbose=False, draw=True, min_size=MIN_SIZE,
                   prefix=PROPHAGE_PREFIX, delimiter=PROPHAGE_DELIMITER,
                   extend_by=EXTEND_BY, run_mode=DEFAULT_RUN_MODE):
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
    :param tmp_dir: path where temporary files can go
    :type tmp_dir: pathlib.Path
    :param cpus: the maximum number of processors to use
    :type cpus: int
    :param verbose: should progress messages be printed along the way?
    :type verbose: bool
    :param draw: should genome diagrams be created at the end?
    :type draw: bool
    :param prefix: prefix for the names of predicted prophage regions
    :type prefix: str
    :param delimiter: delimiter for the prophage region names and numbers
    :type delimiter: str
    :param extend_by: number of basepairs to extend predicted prophage regions
    :type extend_by: int
    :param run_mode: Run mode to operate the prophage identification process
    :type run_mode:
    """
    if verbose:
        print("\tLoading FASTA file...")

    # Parse FASTA file - only keep contigs longer than MIN_LENGTH
    contigs = [x for x in SeqIO.parse(fasta, "fasta") if len(x) >= MIN_LENGTH]

    if verbose:
        print("\tAnnotating CDS and t(m)RNA features...")

    # Set up directory where we can do annotation
    annotation_dir = tmp_dir.joinpath("annotation")
    if not annotation_dir.is_dir():
        annotation_dir.mkdir()

    # Annotate CDS/tRNA/tmRNA features on contigs
    for contig in contigs:
        annotate_contig(contig, annotation_dir) 

    if verbose:
        print("\tMasking conserved bacterial features...")

    mmseqs_dir = tmp_dir.joinpath("mmseqs")
    if not mmseqs_dir.is_dir():
        mmseqs_dir.mkdir()

    # Detect conserved bacterial genes for each contig
    bacterial_masks = assemble_bacterial_mask(contigs, BACTERIAL_REF_FASTA,
                                              BACTERIAL_REF_VALUES, mmseqs_dir)

    if verbose:
        print("\tLooking for high-probability prophage regions...")

    # Predict prophage coordinates for each contig
    prophage_preds = list()
    for contig_i, contig in enumerate(contigs):
        contig_pred = predict_prophage_coords(contig, extend_by,
                                              mask=bacterial_masks[contig_i])

        filtered_contig_pred = []
        for pred in contig_pred: 
            if len(range(*pred)) < min_size:
                continue
            
            filtered_contig_pred.append(pred)

        prophage_preds.append(filtered_contig_pred)

    if all([not any(x) for x in prophage_preds]):
        if verbose:
            print(f"\tNo complete prophages found in {str(fasta)}. PHASTER "
                  f"may be able to find partial (dead) prophages.")
        return

    product_threshold = 0
    if run_mode >= 1:
        if verbose:
            print("\tSearching for phage gene homologs...")

        # Set up directory where we can do hhsearch
        hhsearch_dir = tmp_dir.joinpath("hhsearch")
        if not hhsearch_dir.is_dir():
            hhsearch_dir.mkdir()

        # Search for phage gene remote homologs and
        # annotate the bacterial sequence
        if run_mode >= 2: 
            find_homologs(contigs, prophage_preds, EXTENDED_DB, hhsearch_dir,
                          cpus)
            product_threshold = STRICT_PRODUCT_THRESHOLD
        else:
            find_homologs(contigs, prophage_preds, ESSENTIAL_DB, hhsearch_dir,
                          cpus)
            product_threshold = NORMAL_PRODUCT_THRESHOLD

    prophages = load_initial_prophages(contigs, prophage_preds,
                                       product_threshold=product_threshold,
                                       prefix=prefix, delimiter=delimiter)

    if verbose:
        print("\tIdentifying attL/R...")

    # Set up directory where we can do attL/R detection
    att_dir = tmp_dir.joinpath("att_core")
    if not att_dir.is_dir():
        att_dir.mkdir()

    # Detect attachment sites, where possible, for the predicted prophage
    search_space = extend_by * SEARCH_EXTEND_MULT
    detect_att_sites(prophages, BLASTN_DB, search_space ,att_dir)

    if verbose:
        print("\tGenerating final reports...")

    write_prophage_output(outdir, contigs, prophages, draw)


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def load_initial_prophages(contigs, prophage_predictions,
                           product_threshold=NORMAL_PRODUCT_THRESHOLD,
                           prefix=PROPHAGE_PREFIX,
                           delimiter=PROPHAGE_DELIMITER):
    """Creates Prophage objects from initial prophage prediction coordinates
    and their respective parent SeqRecord objects.

    :param contigs: SeqRecord nucleotide sequence objects
    :type contigs: list[Bio.SeqRecord.SeqRecord]
    :param prophage_predictions: Coordinates for predicted prophages
    :type prophage_predictions: list[list]
    :return: Prophage objects that contain putative sequences and coordinates
    :rtype: list
    """
    prophages = []
    for contig_index, contig in enumerate(contigs):
        # Retrieve the contig seqrecord associated with the coordinates
        contig_predictions = prophage_predictions[contig_index]

        prophage_index = 0
        for prophage_coordinates in contig_predictions:
            # Create a prophage ID from the SeqRecord ID
            prophage_id = "".join([prefix, contig.id,
                                   delimiter, str((prophage_index+1))])
            start = prophage_coordinates[0]
            end = prophage_coordinates[1]

            prophage = Prophage(contig, prophage_id, start=start, end=end)
            prophage.update()

            if len(prophage.products) < product_threshold:
                continue


            products = "\n\t".join(prophage.products)
            prophage_index += 1
            prophages.append(prophage)

    return prophages


def get_reference_map_from_sequence(sequence, sequence_name,
                                    reference_db_path, tmp_dir):
    """Maps sequence BLASTn aligned reference genome IDs to their respective
    alignment result data.

    :param sequence: Query sequence to be aligned to the reference database
    :type sequence: str
    :param sequence_name: Name of the query sequence to be aligned
    :type sequence_name: str
    :param reference_db_path: Path to the database of references to search
    :type reference_db_path: pathlib.Path
    :param tmp_dir: Working directory to place BLASTn inputs and outputs
    :type tmp_dir: pathlib.Path
    :return: A map of aligned reference genome IDs to alignment result data
    """
    sequence_path = tmp_dir.joinpath(".".join([sequence_name, "fasta"]))

    # Write the sequence to a fasta file in the temp directory
    write_fasta([sequence_name], [sequence], sequence_path)

    # Try to retrieve reference results for the sequence to the references
    blast_results = blastn(sequence_path, reference_db_path, tmp_dir)

    reference_map = dict()
    for blast_result in blast_results:
        # Checks to see if the sequence reference ID has already been stored
        results = reference_map.get(blast_result["sseqid"], list())
 
        results.append(blast_result)
        reference_map[blast_result["sseqid"]] = results

    return reference_map


def detect_att_sites(prophages, reference_db_path, search_space,
                     tmp_dir, min_kmer_score=5, sort_key=REF_BLAST_SORT_KEY):
    """Detect attachment sites demarcating predicted prophage regions from
    the bacterial contig.

    :param prophages: Predicted prophages
    :type prophages: list
    :param reference_db_path: Path to the database with reference sequences
    :type reference_db_path: pathlib.Path
    :param search_space: Internal length of the prophage to check for att sites
    :type search_space: int
    :param tmp_dir: Path to place result files.
    :type tmp_dir: pathlib.Path
    :param min_kmer_score: Minimum length threshold of attachment sites.
    :type min_kmer_score: int
    """
    for prophage in prophages:
        # Create prophage working directory within temp directory
        working_dir = tmp_dir.joinpath(prophage.id)
        if not working_dir.is_dir():
            working_dir.mkdir() 

        l_seq = str(prophage.seq[:search_space])
        r_seq = str(prophage.seq[-1*search_space:])

        l_name = f"{prophage.id}_left"
        r_name = f"{prophage.id}_right"

        att_data = find_attachment_site(
                                prophage, l_seq, r_seq,
                                reference_db_path, working_dir, sort_key,
                                k=min_kmer_score, l_name=l_name, r_name=r_name)

        if att_data is not None: 
            prophage.set_coordinates(att_data[0], att_data[1])
            prophage.set_att_len(len(att_data[3]))
        prophage.detect_orientation()

        prophage.update()
        prophage.clean_record()

        prophage.parent_record.features.append(prophage.feature)
        prophage.parent_record.features.sort(key=lambda x: x.location.start)


def write_prophage_output(outdir, contigs, prophages, draw):
    """Generates output structure and writes data to file

    :param outdir: Root directory the data will be written to
    :type outdir: pathlib.Path
    :param contigs: Auto-annotated contigs to be written to file
    :type contigs: list
    :param prophages: Identified prophages to be written to file
    :type prophages: list
    """
    for contig in contigs:
        name = contig.id
        contig.annotations = ANNOTATIONS
        
        genbank_filename = outdir.joinpath(name).with_suffix(".gbk")

        SeqIO.write(contig, genbank_filename, "genbank")
    
    prophage_data_dicts = []
    for prophage in prophages:
        name = prophage.id
        prophage_outdir = outdir.joinpath(name)
        prophage_outdir.mkdir(exist_ok=True)

        genbank_filename = prophage_outdir.joinpath(name).with_suffix(".gbk")
        fasta_filename = prophage_outdir.joinpath(name).with_suffix(".fasta")

        SeqIO.write(prophage.record, genbank_filename, "genbank")
        SeqIO.write(prophage.record, fasta_filename, "fasta")
        if draw:
            diagram_filename = prophage_outdir.joinpath(
                                                    name).with_suffix(".pdf")
            prophage_diagram(prophage.record, diagram_filename)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main(sys.argv[1:])
