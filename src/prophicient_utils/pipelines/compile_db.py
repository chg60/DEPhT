"""Pipeline for downloading, clustering, and making hhsuite3-compatible
databases from phage genomes."""
import argparse
from pathlib import Path

from Prophicient.utilities import (config_handling, entrez, hhsuite,
                                   path_basic)

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = "PhageGeneRef"


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the compile db pipeline

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    args = parse_compile_db(unparsed_args_list)

    config = config_handling.build_complete_config(args.config_file)
    config

    taxon_ids = parse_taxon_ids(args.identifiers_file)

    execute_compile_db(taxon_ids, folder_path=args.folder_path,
                       folder_name=args.folder_name, verbose=args.verbose,
                       force=args.force, cores=args.number_processes)
    pass


def parse_compile_db(unparsed_args_list):
    """Parses export_db arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :typ unparsed_args_list: list[str]
    :returns: Argparse module parsed args.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("identifiers_file", type=path_basic.convert_file_path)

    parser.add_argument("-m", "--folder_name", type=str)
    parser.add_argument("-o", "--folder_path", type=Path)
    parser.add_argument("-c", "--config_file",
                        type=path_basic.convert_file_path)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-f", "--force", action="store_true")
    parser.add_argument("-np", "--number_processes", type=int)

    parser.set_defaults(folder_path=None, folder_name=DEFAULT_FOLDER_NAME,
                        config_file=None, verbose=False, force=False,
                        number_processes=1)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_compile_db(taxon_ids, folder_path=None,
                       folder_name=DEFAULT_FOLDER_NAME, verbose=False,
                       force=False, cores=1):
    working_dir = path_basic.create_working_path(folder_path, folder_name,
                                                 force=force)

    records = retrieve_genome_records(taxon_ids)
    pass


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def parse_taxon_ids(identifiers_file):
    """Parses a taxon identifiers file, assuming each taxon ID is on
    a separate line.

    :param identifiers_file: Path to a file containing taxon identifiers.
    :type identifiers_file: pathlib.Path
    :returns: A set of taxon identifiers.
    :rtype: set()
    """
    taxon_ids = set()

    with identifiers_file.open(mode="r") as filehandle:
        for line in filehandle:
            taxon_ids.add(line.rstrip())

    return taxon_ids


def esearch_taxa(taxon_id, db="nucleotide", idtype="acc", retmax=99999):
    """Searches for genome accession IDs using a taxon identifier.

    :param taxon_id: Taxon identifier to find corresponding accession IDs for.
    :type taxon_id: str
    :param db: Name of the NCBI database to search:
    :type db: str
    :param idtype: Type of identifier to return in the search:
    :type idtype: str
    :param retmax: Number of maximum identifiers to retrieve.
    :type retmax: int
    :returns: Returns specified identifier type related to the taxon ID
    :rtype: list[str]
    """
    taxon_identifier_ref = "".join(["txid", str(taxon_id), "[Orgn]"])

    record = entrez.run_esearch(db=db, idtype=idtype, retmax=retmax,
                                term=taxon_identifier_ref)

    return record["IdList"]


def retrieve_genome_records(taxon_ids):
    """Retrieves genome Biopython SeqRecords from a list of taxon NCBI
    identifiers.

    :param taxon_ids: Taxon identifiers to find corresponding genomes from
    :type taxon_ids: list[str]
    :returns: A list of SeqRecord objects retrieved using taxon identifiers.
    :rtype: list[Bio.SeqRecord.SeqRecord]
    """
    genome_accessions = set()
    for taxon_id in taxon_ids:
        taxon_accessions = esearch_taxa(taxon_id)
        for taxon_accession in taxon_accessions:
            genome_accessions.append(taxon_accession)

    genome_accessions = list(genome_accessions)
    records = entrez.get_records(genome_accessions)

    return records
