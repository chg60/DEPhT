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
from tempfile import mkdtemp

from Bio import SeqIO

from prophicient import PACKAGE_DIR
from prophicient.classes.prophage import ANNOTATIONS, DEFAULT_PRODUCT, Prophage
from prophicient.functions.annotation import annotate_contig, MIN_LENGTH, \
                                             MIN_CDS_FEATURES
from prophicient.functions.att import find_attachment_site
from prophicient.functions.find_homologs import find_homologs
from prophicient.functions.mmseqs import assemble_bacterial_mask
from prophicient.functions.multiprocess import PHYSICAL_CORES
from prophicient.functions.prophage_prediction import predict_prophage_coords
from prophicient.functions.visualization import draw_complete_diagram

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
TMP_DIR = pathlib.Path("/tmp/prophicient")

PROPHAGE_PREFIX = "prophi"
PROPHAGE_DELIMITER = "-"

EXTEND_BY = 4000
ATT_SENSITIVITY = 7
MIN_SIZE = 10000
REF_BLAST_SORT_KEY = "bitscore"

NORMAL_PRODUCT_THRESHOLD = 5
STRICT_PRODUCT_THRESHOLD = 10


# DEFAULT FILE PATHS
# -----------------------------------------------------------------------------
BLASTN_DB = PACKAGE_DIR.joinpath("data/blastn/mycobacteria")

ESSENTIAL_DB = PACKAGE_DIR.joinpath("data/hhsearch/essential/essential")
EXTENDED_DB = PACKAGE_DIR.joinpath("data/hhsearch/extended/extended")

BACTERIAL_REF_FASTA = PACKAGE_DIR.joinpath("data/mmseqs/bacterial_genes.fasta")
BACTERIAL_REF_VALUES = PACKAGE_DIR.joinpath("data/mmseqs/bacterial_genes.pbv")


def parse_args():
    """
    Parse command line arguments.

    :return: parsed_args
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--infile",
                   type=pathlib.Path, nargs="+", required=True,
                   help="path to genome file(s) to scan for prophages")
    p.add_argument("-f", "--input-format",
                   choices=("fasta", "genbank"), default="fasta", type=str,
                   help="input which format your input file(s) are in")
    p.add_argument("-o", "--outdir",
                   type=pathlib.Path, required=True,
                   help="path where outputs can be written")
    p.add_argument("-c", "--cpus",
                   default=PHYSICAL_CORES, type=int,
                   help=f"number of CPU cores to use [default: "
                        f"{PHYSICAL_CORES}]")
    p.add_argument("-n", "--no-draw",
                   action="store_false", dest="draw",
                   help="don't draw genome diagram for identified prophage(s)")
    p.add_argument("-m", "--mode",
                   choices=("fast", "normal", "strict"), default="normal",
                   help="select a runmode that favors speed or accuracy")
    p.add_argument("-s", "--att_sensitivity", type=float,
                   default=ATT_SENSITIVITY,
                   help="sensitivity parameter for att site detection.")
    p.add_argument("-d", "--dump-data", action="store_true",
                   help="dump all support data to outdir")
    p.add_argument("-v", "--verbose",
                   action="store_true",
                   help="print progress messages as the program runs")

    return p.parse_args()


def main():
    """
    Main function that interfaces with command line args and the
    program workflow.
    """
    args = parse_args()

    # Get input file(s) and format
    infiles, fmt = args.infile, args.input_format

    # Are we going to draw genome diagrams?
    draw = args.draw

    # Are we dumping data into outdir when done?
    dump = args.dump_data

    # What runmode are we using - never draw in fast runmode
    runmode = args.mode
    if runmode == "fast" and draw:
        # Why do we want to turn off draw?
        # print("fast mode - turning off draw...")
        # draw = False
        pass

    # Print verbose outputs?
    verbose = args.verbose

    # How many physical cores can we use?
    cpus = args.cpus

    # What's the att sensitivity modifier?
    att_sensitivity = args.att_sensitivity

    # Get output dir and make sure it's a valid path
    outdir = pathlib.Path(args.outdir).resolve()
    if not outdir.is_dir():
        print(f"'{str(outdir)}' does not exist - creating it...")
        outdir.mkdir(parents=True)

    # Refresh the temporary directory
    if TMP_DIR.is_dir():
        shutil.rmtree(TMP_DIR)
    TMP_DIR.mkdir(parents=True)

    # Mark program start time
    mark = datetime.now()

    # OK, let's go!
    for infile in infiles:
        # Make sure filepath exists - skip if it doesn't
        if not infile.is_file():
            print(f"'{str(infile)}' does not exist - skipping it...")
            continue

        # Skip .DS_Store files
        if infile.name == ".DS_Store":
            print("skipping .DS_Store file...")
            continue

        # Set up a temporary directory for this genome
        tmp_dir = pathlib.Path(mkdtemp(dir=TMP_DIR))
        if not tmp_dir.is_dir():
            tmp_dir.mkdir()

        # Parse all contigs of annotation-worthy length
        if verbose:
            print(f"\nparsing '{str(infile)}'...")
        contigs = [x for x in SeqIO.parse(infile, fmt) if len(x) >= MIN_LENGTH]

        if not contigs:
            print(f"no {fmt}-formatted records found in '{str(infile)}' - "
                  f"skipping it...")
            shutil.rmtree(tmp_dir)  # clean up after ourselves
            continue

        # Annotate contigs if format is "fasta"
        if fmt == "fasta":
            if verbose:
                print("annotating t(m)RNA and CDS genes de novo...")

            annotate_dir = tmp_dir.joinpath("annotate")
            if not annotate_dir.is_dir():
                annotate_dir.mkdir()

            for contig in contigs:
                annotate_contig(contig, annotate_dir)

        else:
            if verbose:
                print("using flat file annotation...")

            load_flat_file_contigs(contigs)

        # Filter contigs that don't have enough CDS features
        contigs = [contig for contig in contigs
                   if (len([x for x in contig.features if x.type == "CDS"]) >=
                       MIN_CDS_FEATURES)]

        if not contigs:
            print(f"no contigs long enough to analyze in '{str(infile)}' - "
                  f"skipping it...")
            shutil.rmtree(tmp_dir)  # clean up after ourselves
            continue

        if verbose:
            print("masking conserved bacterial features...")

        mmseqs_dir = tmp_dir.joinpath("mmseqs")
        if not mmseqs_dir.is_dir():
            mmseqs_dir.mkdir()

        # Detect conserved bacterial genes for each contig
        bacterial_masks = assemble_bacterial_mask(
            contigs, BACTERIAL_REF_FASTA, BACTERIAL_REF_VALUES, mmseqs_dir)

        if verbose:
            print("looking for high-probability prophage regions...")

        # Predict prophage coordinates for each contig
        prophage_preds = list()
        for i, contig in enumerate(contigs):
            contig_pred = predict_prophage_coords(contig, EXTEND_BY,
                                                  mask=bacterial_masks[i])

            filtered_contig_pred = []
            for pred in contig_pred:
                if len(range(*pred)) < MIN_SIZE:
                    continue

                filtered_contig_pred.append(pred)

            prophage_preds.append(filtered_contig_pred)

        if all([not any(x) for x in prophage_preds]):
            print(f"no complete prophages found in {str(infile)}. "
                  f"PHASTER may be able to find partial (dead) prophages.")
            shutil.rmtree(tmp_dir)  # clean up after ourselves
            continue

        product_threshold = 0
        if runmode in ("normal", "strict"):
            if verbose:
                print("searching for phage gene homologs...")

            hhsearch_dir = tmp_dir.joinpath("hhsearch")
            if not hhsearch_dir.is_dir():
                hhsearch_dir.mkdir()

            # Search for phage gene remote homologs
            if runmode == "strict":
                find_homologs(contigs, prophage_preds, EXTENDED_DB,
                              hhsearch_dir, cpus)
                product_threshold = STRICT_PRODUCT_THRESHOLD
            else:
                find_homologs(contigs, prophage_preds, ESSENTIAL_DB,
                              hhsearch_dir, cpus)
                product_threshold = NORMAL_PRODUCT_THRESHOLD

        prophages = load_initial_prophages(contigs, prophage_preds,
                                           product_threshold,
                                           prefix=PROPHAGE_PREFIX,
                                           delimiter=PROPHAGE_DELIMITER)

        if not prophages:
            print(f"no complete prophages found in {str(infile)}. "
                  f"PHASTER may be able to find partial (dead) prophages.")
            shutil.rmtree(tmp_dir)  # clean up after ourselves
            continue

        if verbose:
            print("searching for attL/R...")

        # Set up directory where we can do attL/R detection
        att_dir = tmp_dir.joinpath("att_core")
        if not att_dir.is_dir():
            att_dir.mkdir()

        # Detect attachment sites, where possible, for the predicted prophage
        search_space = int(EXTEND_BY * att_sensitivity)
        detect_att_sites(prophages, BLASTN_DB, search_space, att_dir)

        if verbose:
            print("generating final reports...")

        genome_outdir = outdir.joinpath(f"{infile.stem}")
        if not genome_outdir.is_dir():
            genome_outdir.mkdir()

        draw_dir = tmp_dir.joinpath("draw_diagram")
        if not draw_dir.is_dir():
            draw_dir.mkdir()
        write_prophage_output(genome_outdir, contigs, prophages, draw_dir,
                              draw)

        if dump:
            destination = genome_outdir.joinpath("tmp_data")
            shutil.copytree(tmp_dir, destination)
        shutil.rmtree(tmp_dir)  # clean up after ourselves

    print(f"\nTotal runtime: {str(datetime.now() - mark)}")


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def load_flat_file_contigs(contigs):
    """Function to clean up and format SeqRecord sequence contigs created
    from imported flat file annotations.

    :param contigs: Imported sequence contigs
    :type contigs: list
    """
    for contig in contigs:
        for feature in contig.features:
            if feature.type != "CDS":
                continue

            if not feature.qualifiers.get("translation"):
                dna = feature.extract(contig.seq)
                translation = dna.translate(to_stop=True, table=11)
                feature.qualifiers["translation"] = [str(translation)]

            feature.qualifiers["product"] = [DEFAULT_PRODUCT]


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
            if start < 0:
                start = 0

            end = prophage_coordinates[1]
            if end >= len(contig.seq):
                end = len(contig.seq) - 1

            prophage = Prophage(contig, prophage_id, start=start, end=end)
            prophage.update()

            if len(prophage.products) < product_threshold:
                continue

            prophage_index += 1
            prophages.append(prophage)

    return prophages


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

        half_len = len(prophage.seq) // 2
        if search_space > half_len:
            search_space = half_len

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


def write_prophage_output(outdir, contigs, prophages, tmp_dir, draw):
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

    for prophage in prophages:
        name = prophage.id
        prophage_outdir = outdir.joinpath(name)
        prophage_outdir.mkdir(exist_ok=True)

        genbank_filename = prophage_outdir.joinpath(name).with_suffix(".gbk")
        fasta_filename = prophage_outdir.joinpath(name).with_suffix(".fasta")

        SeqIO.write(prophage.record, genbank_filename, "genbank")
        SeqIO.write(prophage.record, fasta_filename, "fasta")

    if draw:
        draw_complete_diagram(outdir, contigs, prophages, tmp_dir,
                              name=outdir.stem)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    main()
