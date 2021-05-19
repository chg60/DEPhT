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
from prophicient.classes.prophage import Prophage
from prophicient.functions.annotation import annotate_contig, MIN_LENGTH
from prophicient.functions.att import find_attachment_site
from prophicient.functions.blastn import blastn
from prophicient.functions.fasta import write_fasta
from prophicient.functions.find_homologs import find_homologs
from prophicient.functions.multiprocess import PHYSICAL_CORES
from prophicient.functions.prophage_prediction import predict_prophage_coords
from prophicient.functions.visualization import prophage_diagram

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
RUN_MODE_MAP = {"fast": 0, "norm": 1, "full": 2}
DEFAULT_RUN_MODE = "norm"

TMP_DIR = pathlib.Path("/tmp/prophicient")

BLASTN_DB = PACKAGE_DIR.joinpath("data/blastn/mycobacteria")
HHSEARCH_DB = PACKAGE_DIR.joinpath("data/hhsearch/functions")

EXTEND_BY = 5000
REF_BLAST_SORT_KEY = "bitscore"
PRODUCT_THRESHOLD = 3


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
    parser.add_argument("--verbose", action="store_true",
                        help="toggles verbosity of pipeline")
    parser.add_argument("--no-draw", action="store_true",
                        help="don't output genome map PDFs for identified "
                             "prophages")
    parser.add_argument("--mode", type=str,
                        choices=RUN_MODE_MAP.keys(), default=DEFAULT_RUN_MODE,
                        help="Choose run mode of the pipeline.")
    parser.add_argument("--cpus", type=int, default=PHYSICAL_CORES,
                        help=f"number of processors to use [default: "
                             f"{PHYSICAL_CORES}]")
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

    # Create temporary dir, if it doesn't exist.
    if not TMP_DIR.is_dir():
        TMP_DIR.mkdir()

    cpus = args.cpus
    verbose = args.verbose
    draw = not args.no_draw
    mode = RUN_MODE_MAP[args.mode]

    # Mark program start time
    mark = datetime.now()
    find_prophages(infile, outdir, TMP_DIR, cpus, verbose, draw, EXTEND_BY,
                   mode)
    print(f"\nTotal runtime: {str(datetime.now() - mark)}")

    # Clean up after ourselves
    shutil.rmtree(TMP_DIR)


def find_prophages(fasta, outdir, tmp_dir, cpus, verbose, draw, extend_by,
                   run_mode):
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
        print("\tLooking for high-probability prophage regions...")

    # Predict prophage coordinates for each contig
    prophage_preds = list()
    for contig in contigs:
        prophage_preds.append(predict_prophage_coords(contig, EXTEND_BY))

    if all([not any(x) for x in prophage_preds]):
        if verbose:
            print(f"\tNo complete prophages found in {str(fasta)}. PHASTER "
                  f"may be able to find partial (dead) prophages.")
        return

    if run_mode >= 1:
        if verbose:
            print("\tSearching for phage gene homologs...")

        # Set up directory where we can do hhsearch
        hhsearch_dir = tmp_dir.joinpath("hhsearch")
        if not hhsearch_dir.is_dir():
            hhsearch_dir.mkdir()

        # Search for phage gene remote homologs and
        # annotate the bacterial sequence
        find_homologs(contigs, prophage_preds, HHSEARCH_DB, hhsearch_dir, cpus)

    prophages = load_initial_prophages(contigs, prophage_preds,
                                       reject=(run_mode >=1))

    if verbose:
        print("\tSearching for attL/R...")

    # Set up directory where we can do attL/R detection
    att_dir = tmp_dir.joinpath("att_core")
    if not att_dir.is_dir():
        att_dir.mkdir()

    # Detect attachment sites, where possible, for the predicted prophage
    detect_att_sites(prophages, BLASTN_DB, extend_by*2, att_dir)

    if verbose:
        print("\tGenerating final reports...")

    write_prophage_output(outdir, contigs, prophages, draw)


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def load_initial_prophages(contigs, prophage_predictions, reject=True,
                           product_threshold=PRODUCT_THRESHOLD):
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
            prophage_id = "".join(["prophi", contig.id,
                                   "-", str((prophage_index+1))])
            start = prophage_coordinates[0]
            end = prophage_coordinates[1]

            prophage = Prophage(contig, prophage_id, start=start, end=end)
            prophage.update()

            if reject:
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
        result_exists = reference_map.get(blast_result["sseqid"], False)

        if not result_exists:
            reference_map[blast_result["sseqid"]] = blast_result

    return reference_map


def detect_att_sites(prophages, reference_db_path, extend_by,
                     tmp_dir, min_kmer_score=5, sort_key=REF_BLAST_SORT_KEY):
    """Detect attachment sites demarcating predicted prophage regions from
    the bacterial contig.

    :param prophages: Predicted prophages
    :type prophages: list
    :param reference_db_path: Path to the database with reference sequences
    :type reference_db_path: pathlib.Path
    :param extend_by: Internal length of the prophage to check for att sites
    :type extend_by: int
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

        # If the internal extend_by will cause the right and left regions
        # to overlap, limit the external extend_by to half the region
        seq_len = len(prophage.seq)
        half_len = seq_len // 2
        if extend_by > half_len:
            extend_by = half_len

        # BLASTn the left region against the reference database
        left_seq = str(prophage.seq[:half_len])
        left_name = f"{prophage.id}_left"
        left_map = get_reference_map_from_sequence(
            left_seq, left_name, reference_db_path, working_dir)

        # BLASTn the right region against the reference database
        right_seq = str(prophage.seq[half_len:])
        right_name = f"{prophage.id}_right"
        right_map = get_reference_map_from_sequence(
            right_seq, right_name, reference_db_path, working_dir)

        # Find BLASTn reference genomes that both regions aligned to
        ref_ids = list(set(left_map.keys()).intersection(right_map.keys()))

        # Sort reference genomes by best cumulative E-value
        ref_data = \
            [(left_map[ref_id], right_map[ref_id]) for ref_id in ref_ids]
        ref_data.sort(key=lambda x:
                            (float(x[0][sort_key]) + float(x[1][sort_key])),
                      reverse=True)

        new_coords = None
        for data_tuple in ref_data:
            l_data = data_tuple[0]
            r_data = data_tuple[1]

            # Find the coordinate ranges of the aligned reference genome
            left_ref_range = range(int(l_data["sstart"]),
                                   int(l_data["send"]))
            right_ref_range = range(int(r_data["sstart"]),
                                    int(r_data["send"]))

            # If the coordinates indicate overlap, determine the overlap length
            overlap_range = set(left_ref_range).intersection(
                                                    set(right_ref_range))
            overlap_len = len(overlap_range)

            # If the overlap length meets the minimum att length
            # treat the sequence as a putative attB that we can use 
            # as a reference to set the boundaries for the predicted prophage
            if overlap_len >= min_kmer_score:
                # Find the right coordinate of the putative attL
                # in the aligned left region of the prophage 
                l_qend = int(l_data["qend"])

                # Find the left coordinate of the putative attR
                # in the aligned right region of the prophage 
                r_qstart = int(r_data["qstart"]) 

                # Set the prophage start as the left coordinate of the attL
                # of the putative prophage added to the start of the left
                # region of the prophage in the bacterial sequence
                new_start = (prophage.start + int(l_qend) - overlap_len)
                # Set the prophage end as the right coordinate of the attR
                # of the putative prophage added to the start of the right
                # region of the prophage in the bacterial sequence
                new_end = ((prophage.end - half_len) +
                            int(r_qstart) + overlap_len)
                
                new_coords = (new_start, new_end)
                att_len = overlap_len
                break

        if not new_coords:
            # Sets the putative origin to half the internal extensions
            l_origin = extend_by // 2
            r_origin = extend_by // 2
            # If any part of the region aligned, use the last aligned
            # coordinate as the origin
            if ref_data:
                data_tuple = ref_data[0]

                l_end = int(data_tuple[0]["qend"])
                if l_end < l_origin:
                    l_origin = l_end

                r_start = int(data_tuple[1]["qstart"])
                if r_start > r_origin:
                    r_origin = r_start
            
            # Attempts to find an attachment site by finding the longest
            # sequence represented in both the left and right prophage regions,
            # penalizing length with distance from the set origin coordinate
            kmer_data = find_attachment_site(prophage, left_seq, right_seq,
                                             l_origin, r_origin, working_dir,
                                             l_name=left_name,
                                             r_name=right_name,
                                             k=min_kmer_score)

            new_coords = (prophage.start, prophage.end)
            att_len = 0
            if kmer_data is not None:
                if kmer_data[2] >= min_kmer_score:
                    new_start = (prophage.start +
                                 kmer_data[0].location.start)
                    new_end = ((prophage.end - half_len) +
                               kmer_data[1].location.end)
                    new_coords = (new_start, new_end)
                    att_len = len(kmer_data[3])

        prophage.set_coordinates(*new_coords)
        prophage.set_att_len(att_len)
        prophage.detect_orientation()

        prophage.update()
        prophage.clean_record()


def write_prophage_output(outdir, contigs, prophages, draw):
    """Generates output structure and writes data to file

    :param outdir: Root directory the data will be written to
    :type outdir: pathlib.Path
    :param contigs: Auto-annotated contigs to be written to file
    :type contigs: list
    :param prophages: Identified prophages to be written to file
    :type prophages: list
    """
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
