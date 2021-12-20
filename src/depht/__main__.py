"""
DEPhT scans bacterial genomes looking for prophages. Regions
identified as prophage candidates are further scrutinized, and
attachment sites identified as accurately as possible before
prophage extraction and generating the final report.
"""
import argparse
import csv
import pathlib
import shutil
import sys
from datetime import datetime
from tempfile import mkdtemp

from Bio import SeqIO

from depht.classes.contig import CODING_FEATURE_TYPES, Contig
from depht.classes.prophage import ANNOTATIONS, DEFAULT_PRODUCT, Prophage
from depht.functions.annotation import annotate_record, MIN_LENGTH
from depht.functions.att import find_attachment_site
from depht.functions.find_homologs import find_homologs
from depht.functions.mmseqs import assemble_bacterial_mask
from depht.functions.multiprocess import PHYSICAL_CORES
from depht.functions.prophage_prediction import predict_prophage_coords, WINDOW
from depht.functions.sniff_format import sniff_format
from depht.functions.visualization import draw_complete_diagram

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEPHT_DIR = pathlib.Path().home().joinpath(".depht")
if not DEPHT_DIR.is_dir():
    DEPHT_DIR.mkdir()

MODEL_DIR = DEPHT_DIR.joinpath("models")
if not MODEL_DIR.is_dir():
    MODEL_DIR.mkdir()

LOCAL_MODELS = [x.name for x in MODEL_DIR.iterdir() if x.name != ".DS_Store"]
if len(LOCAL_MODELS) == 0:
    # TODO: automatically retrieve models from OSF if the user says we can
    print("No DEPhT models found in ~/.depht/models. Please create one using "
          "'depht_train', or download models from https://osf.io/zt4n3/")
    sys.exit(0)

# Where temporary file I/O should happen
TMP_DIR = DEPHT_DIR.joinpath("tmp")

# Model can't scan contigs with fewer CDS than window size
MIN_CDS_FEATURES = WINDOW

# For naming any identified prophages
PROPHAGE_PREFIX = "prophi"
PROPHAGE_DELIMITER = "-"

# For attachment site detection
EXTEND_BY = 5000
ATT_SENSITIVITY = 7
REF_BLAST_SORT_KEY = "bitscore"

# For deciding whether to cull predicted "prophages"
MIN_SIZE = 20000
MIN_PRODUCTS_NORMAL = 5
MIN_PRODUCTS_STRICT = 10

CONTIG_DATA_HEADER = ["Gene ID", "Start", "End", "Prediction",
                      "Bacterial Homology", "Phage Homology"]


def parse_args():
    """
    Parse command line arguments.

    :return: parsed_args
    """
    p = argparse.ArgumentParser(description=__doc__, prog="depht")
    # Positional arguments
    p.add_argument("infile", type=pathlib.Path, nargs="+",
                   help="path to genome file(s) to scan for prophages")
    p.add_argument("outdir",  type=pathlib.Path,
                   help="path where output files can be written")

    # Optional arguments
    p.add_argument("--model", type=str,
                   choices=set(LOCAL_MODELS), default=LOCAL_MODELS[0],
                   help=f"which local model should be used [default: "
                        f"{LOCAL_MODELS[0]}]")
    p.add_argument("-c", "--cpus",
                   metavar="", default=PHYSICAL_CORES, type=int,
                   help=f"number of CPU cores to use [default: "
                        f"{PHYSICAL_CORES}]")
    p.add_argument("-n", "--no-draw", action="store_false", dest="draw",
                   help="don't draw genome diagram for identified prophage(s)")
    p.add_argument("-m", "--mode",
                   choices=("fast", "normal", "sensitive"), default="normal",
                   help="select a runmode that favors speed or accuracy")
    p.add_argument("-s", "--att_sens", type=float,
                   default=ATT_SENSITIVITY, metavar="",
                   help="sensitivity parameter for att site detection.")
    p.add_argument("-d", "--dump-data", action="store_true",
                   help="dump all support data to outdir")
    p.add_argument("-v", "--verbose",
                   action="store_true",
                   help="print progress messages as the program runs")
    p.add_argument("-t", "--tmp-dir",
                   type=pathlib.Path, default=TMP_DIR, metavar="",
                   help=f"temporary directory to use for file I/O "
                        f"[default: {TMP_DIR}]")
    p.add_argument("-p", "--products",
                   type=int, default=None, metavar="",
                   help="minimum number of phage homologs to report a prophage")
    p.add_argument("-l", "--length",
                   type=int, default=MIN_SIZE, metavar="",
                   help=f"minimum length to report a prophage "
                        f"[default: {MIN_SIZE}]")

    return p.parse_args()


def main():
    """
    Main function that interfaces with command line args and the
    program workflow.
    """
    # Invoke the help menu if no args given
    if len(sys.argv) == 1:
        sys.argv.append("-h")

    args = parse_args()

    infiles = args.infile       # What input file(s) did the user give?
    model = args.model          # Which model was selected?
    draw = args.draw            # Are we going to draw genome diagrams?
    dump = args.dump_data       # Are we dumping data into outdir when done?
    runmode = args.mode         # What runmode are we using?
    verbose = args.verbose      # Print verbose outputs?
    cpus = args.cpus            # How many physical cores can we use?
    att_sens = args.att_sens    # What's the att sensitivity modifier?
    min_length = args.length    # Minimum prophage length in bp

    # Get output dir and make sure it's a valid path
    outdir = pathlib.Path(args.outdir).resolve()
    if not outdir.is_dir():
        print(f"'{str(outdir)}' does not exist - creating it...")
        outdir.mkdir(parents=True)

    # Get the temporary directory, refresh it if it exists
    tmpdir = pathlib.Path(args.tmp_dir).resolve()
    if not tmpdir.is_dir():
        tmpdir.mkdir(parents=True)

    # Set up model directories
    model_dir = MODEL_DIR.joinpath(model)
    bact_ref_fasta = model_dir.joinpath("mmseqs/bacterial_genes.fasta")
    bact_ref_values = model_dir.joinpath("mmseqs/bacterial_genes.pbv")
    essential_db = model_dir.joinpath("hhsearch/essential")
    extended_db = model_dir.joinpath("hhsearch/extended")
    blast_db = model_dir.joinpath("blastn/references")
    classifier = model_dir.joinpath("classifier.pkl")

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
        genome_tmp_dir = pathlib.Path(mkdtemp(dir=tmpdir))
        if not genome_tmp_dir.is_dir():
            genome_tmp_dir.mkdir()

        # Automatically check the file format
        fmt = sniff_format(infile)

        # Parse all contigs of annotation-worthy length
        if verbose:
            print(f"\nparsing '{str(infile)}' as {fmt}...")

        # Parse the file and retain contigs above cutoff length
        records = [x for x in SeqIO.parse(infile, fmt) if len(x) >= MIN_LENGTH]

        if not records:
            if verbose:
                print(f"no {fmt}-formatted records of at least {MIN_LENGTH}bp "
                      f"found in '{str(infile)}' - skipping it...")

            shutil.rmtree(genome_tmp_dir)  # clean up after ourselves
            continue

        # Annotate contigs if format is "fasta"
        if fmt == "fasta":
            if verbose:
                print("annotating t(m)RNA and CDS genes de novo...")

            annotate_dir = genome_tmp_dir.joinpath("annotate")
            if not annotate_dir.is_dir():
                annotate_dir.mkdir()

            for record in records:
                annotate_record(record, annotate_dir)

        else:
            if verbose:
                print("using flat file annotation...")

            cleanup_flatfile_records(records)

        # Filter contigs that don't have enough CDS features
        records = [record for record in records if (len(
            [x for x in record.features if x.type == "CDS"]) >=
                                                    MIN_CDS_FEATURES)]

        if not records:
            print(f"no contigs with enough CDS features to analyze in "
                  f"'{str(infile)}' - skipping it...")

            shutil.rmtree(genome_tmp_dir)  # clean up after ourselves
            continue

        contigs = load_contigs(records)

        if verbose:
            print("masking conserved bacterial features...")

        mmseqs_dir = genome_tmp_dir.joinpath("mmseqs")
        if not mmseqs_dir.is_dir():
            mmseqs_dir.mkdir()

        # Detect conserved bacterial genes for each contig
        bacterial_masks = assemble_bacterial_mask(
            contigs, bact_ref_fasta, bact_ref_values, mmseqs_dir)

        if verbose:
            print("looking for high-probability prophage regions...")

        # Predict prophage coordinates for each contig
        prophage_predictions = list()
        for i, contig in enumerate(contigs):
            prediction = predict_prophage_coords(
                contig, classifier, EXTEND_BY, mask=bacterial_masks[i])

            filtered_prediction = []
            for pred in prediction:
                if len(range(*pred)) < min_length:
                    continue

                filtered_prediction.append(pred)

            prophage_predictions.append(filtered_prediction)

        if all([not any(x) for x in prophage_predictions]) and not dump:
            print(f"no complete prophages found in '{str(infile)}'...")

            shutil.rmtree(genome_tmp_dir)  # clean up after ourselves
            continue

        product_threshold = 0
        if runmode in ("normal", "sensitive"):
            hhsearch_dir = genome_tmp_dir.joinpath("hhsearch")
            if not hhsearch_dir.is_dir():
                hhsearch_dir.mkdir()

            # Search for phage gene remote homologs
            if verbose:
                print("searching for phage gene homologs...")

            find_homologs(contigs, prophage_predictions, essential_db,
                          hhsearch_dir, cpus)
            product_threshold = MIN_PRODUCTS_NORMAL

            if runmode == "sensitive":
                if verbose:
                    print("extending search for phage gene homologs...")
                find_homologs(contigs, prophage_predictions, extended_db,
                              hhsearch_dir, cpus, cache_scores=False)
                product_threshold = MIN_PRODUCTS_STRICT
        else:
            for contig in contigs:
                contig.fill_hhsearch_scores()

        if args.products is not None:
            product_threshold = args.products

        prophages = load_initial_prophages(contigs, prophage_predictions,
                                           product_threshold,
                                           prefix=PROPHAGE_PREFIX,
                                           delimiter=PROPHAGE_DELIMITER)

        if verbose and prophages:
            print("searching for attL/R...")

        # Set up directory where we can do attL/R detection
        att_dir = genome_tmp_dir.joinpath("att_core")
        if not att_dir.is_dir():
            att_dir.mkdir()

        # Detect attachment sites, where possible, for the predicted prophage
        search_space = att_sens * EXTEND_BY
        detect_att_sites(prophages, blast_db, search_space, att_dir)

        prophages = [prophage for prophage in prophages
                     if prophage.length >= min_length]

        if not prophages and not dump:
            print(f"no complete prophages found in '{str(infile)}'...")

            shutil.rmtree(genome_tmp_dir)  # clean up after ourselves
            continue

        if verbose:
            print("generating final reports...")

        genome_outdir = outdir.joinpath(f"{infile.stem}")
        if not genome_outdir.is_dir():
            genome_outdir.mkdir()

        draw_dir = genome_tmp_dir.joinpath("draw_diagram")
        if not draw_dir.is_dir():
            draw_dir.mkdir()
        write_prophage_output(genome_outdir, contigs, prophages, draw_dir,
                              draw)

        if dump:
            destination = genome_outdir.joinpath("tmp_data")

            if destination.exists():
                shutil.rmtree(destination)

            shutil.copytree(genome_tmp_dir, destination)
        shutil.rmtree(genome_tmp_dir)  # clean up after ourselves

    print(f"\nTotal runtime: {str(datetime.now() - mark)}")


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def cleanup_flatfile_records(records):
    """
    Function to clean up and format SeqRecord sequence contigs created
    from imported flat file annotations.

    :param records: imported records
    :type records: list
    """
    for record in records:
        features = list()
        for feature in record.features:
            if feature.type not in CODING_FEATURE_TYPES:
                continue

            if feature.type == "CDS":
                if not feature.qualifiers.get("translation"):
                    dna = feature.extract(record.seq)
                    translation = dna.translate(to_stop=True, table=11)
                    feature.qualifiers["translation"] = [str(translation)]

                feature.qualifiers["product"] = [DEFAULT_PRODUCT]

            features.append(feature)

        record.features = features


def load_contigs(contig_records):
    """Function to create Contig objects from bacterial sequence contig
    SeqRecords

    :param contig_records: Bacterial sequence contig SeqRecords
    :type contig_records: list
    :return: Contig options built from inputted SeqRecords
    :rtype: list
    """
    contigs = list()
    for record in contig_records:
        contig = Contig(record, record.name)
        contigs.append(contig)

    return contigs


def load_initial_prophages(contigs, prophage_predictions,
                           product_threshold=MIN_PRODUCTS_NORMAL,
                           prefix=PROPHAGE_PREFIX,
                           delimiter=PROPHAGE_DELIMITER):
    """Creates Prophage objects from initial prophage prediction coordinates
    and their respective parent SeqRecord objects.

    :param contigs: SeqRecord nucleotide sequence objects
    :type contigs: list[Contig]
    :param prophage_predictions: coordinates for predicted prophages
    :type prophage_predictions: list[list]
    :param product_threshold: number of products
    :type product_threshold:
    :param prefix: how should locus tags begin?
    :type prefix: str
    :param delimiter: how should locus tags be delimited
    :type delimiter: str
    :return: prophage objects that contain putative sequences and coordinates
    :rtype: list
    """
    prophages = []
    for contig_index, contig in enumerate(contigs):
        # Retrieve the contig SeqRecord associated with the coordinates
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

            prophage = Prophage(contig.record, prophage_id, start=start,
                                end=end)
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
    :param sort_key: the blast column to sort possible att pairs on
    :type sort_key: str
    """
    for prophage in prophages:
        # Create prophage working directory within temp directory
        working_dir = tmp_dir.joinpath(prophage.id)
        if not working_dir.is_dir():
            working_dir.mkdir()

        half_len = len(prophage.seq) // 2
        if search_space > half_len:
            search_space = half_len

        search_space = int(search_space)

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
    :param tmp_dir: where this genome's temporary data are found
    :type tmp_dir: pathlib.Path
    :param draw: draw diagram(s) for this genome's prophage(s)?
    :type draw: bool
    """
    for contig in contigs:
        name = contig.id
        contig.record.annotations = ANNOTATIONS

        genbank_filename = outdir.joinpath(f"{name}.gbk")
        table_filename = outdir.joinpath(f"{name}.csv")

        SeqIO.write(contig.record, genbank_filename, "genbank")

        write_contig_data(contig, table_filename)

    for prophage in prophages:
        name = prophage.id
        prophage_outdir = outdir.joinpath(name)
        prophage_outdir.mkdir(exist_ok=True)

        genbank_filename = prophage_outdir.joinpath(f"{name}.gbk")
        fasta_filename = prophage_outdir.joinpath(f"{name}.fasta")

        SeqIO.write(prophage.record, genbank_filename, "genbank")
        SeqIO.write(prophage.record, fasta_filename, "fasta")

    if draw and prophages:
        draw_complete_diagram(outdir, [contig.record for contig in contigs],
                              prophages, tmp_dir, name=outdir.name)


def write_contig_data(contig, outpath):
    """Generates a csv from data associated with each gene from a bacterial
    sequence contig.

    :param contig: Bacterial sequence contig class
    :type contig: prophicient.classes.contig.Contig
    :param outpath: Path to the outputted data table file
    :type outpath: pathlib.Path or str
    """
    handle = open(outpath, "w")
    csv_writer = csv.DictWriter(handle, fieldnames=CONTIG_DATA_HEADER)

    csv_writer.writeheader()

    for i, feature in enumerate(contig.genes):
        data = (contig.gene_ids[i], feature.location.start,
                feature.location.end, contig.model_scores[i],
                contig.mask_bits[i], contig.hhsearch_scores[i])

        data_dict = dict()
        for j, header in enumerate(CONTIG_DATA_HEADER):
            data_dict[header] = data[j]

        csv_writer.writerow(data_dict)

    handle.close()


if __name__ == "__main__":
    main()
