"""
Prophicient scans bacterial genomes looking for prophages. Regions
identified as prophage candidates are further scrutinized, and
attachment sites identified as accurately as possible before
prophage extraction and generating the final report.
"""
import re
import sys
import argparse
import pathlib
import time

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from prophicient.functions.fasta import write_fasta
from prophicient.functions.multiprocess import CPUS
from prophicient.functions.wrapper_basic import autoannotate
from prophicient.functions.prefilter import prefilter_genome
from prophicient.functions.find_homologs import find_homologs

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
PRODIGAL_FORMAT = re.compile(r">\w+_(\d+) # (\d+) # (\d+) # (-?1) # ID=\w+;"
                             r"partial=(\d+);start_type=(\w+);"
                             r"rbs_motif=(.*|None);"
                             r"rbs_spacer=(.*bp|None);gc_cont=(\d.\d+)\s+")

ANNOTATIONS = {"molecule_type": "DNA",
               "topology": "linear",
               "data_file_division": "PHG",
               "date": "",
               "accessions": [],
               "sequence_version": "1",
               "keywords": [],
               "source": "",
               "organism": "",
               "taxonomy": [],
               "comment": [""]}


def parse_args(arguments):
    """
    Parse command line arguments
    :param arguments:  command line arguments that program was invoked with
    :type arguments: list
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--infile", type=pathlib.Path, required=True,
                        help=("FASTA file containing nucleotide sequence to "
                              "scan for prophages"))
    parser.add_argument("-o", "--outdir", type=pathlib.Path,
                        required=True,
                        help="path where output files can be written")
    parser.add_argument("-d", "--database", type=pathlib.Path,
                        required=True,
                        help="Path to a Prophicient database.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Toggles verbosity of pipeline.")
    parser.add_argument("--cpus", type=int, default=CPUS,
                        help=f"number of processors to use [default: {CPUS}]")
    return parser.parse_args(arguments)


def main(arguments):
    args = parse_args(arguments)

    # Verify that the input filepath is valid
    infile = args.infile
    if not infile.is_file():
        print(f"'{str(infile)}' is not a valid input file - exiting")
        sys.exit(1)

    if not args.outdir.is_dir():
        print(f"'{str(args.outdir)}' does not exist - creating it")
        args.outdir.mkdir(parents=True)

    start = time.time()
    execute_prophicient(args.infile, args.outdir, args.database,
                        cores=args.cpus, verbose=args.verbose)
    stop = time.time()

    print("\nTime elapsed: {:.3f}".format(stop - start))


def execute_prophicient(infile, outdir, database, cores=1, verbose=False):
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
            for region_data, integrase_homologs in putative_prophage_regions:
                putative_prophage = extract_prophage(region_data[0],
                                                     integrase_homologs)

                extracted_prophages.append(putative_prophage)


def extract_prophage(record, integrase_features):
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

        integrase_features = find_integrase_homologs(region_tuple[0], data_dir,
                                                     database, cores=cores,
                                                     verbose=verbose)

        if not integrase_features:
            continue

        if verbose:
            for feature in integrase_features:
                print("...Integrase homolog found at "
                      f"{(region_start + feature.location.start)}")

        putative_prophage_regions.append((region_tuple, integrase_features))

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

    phagey_gene_ids = find_homologs(trans_dir, hhresults_dir, database,
                                    cores=cores, verbose=verbose)

    integrase_features = [gene_id_feature_map[gene_id]
                          for gene_id in phagey_gene_ids]

    return integrase_features


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
