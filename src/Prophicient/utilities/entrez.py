from Bio import (Entrez, SeqIO)


def set_entrez_credentials(tool=None, email=None, api_key=None):
    """Set BioPython Entrez credentials to improve speed and reliability.

    :param tool: Name of the software/tool being used.
    :type tool: str
    :param email: Email contact information for NCBI.
    :type email: str
    :param api_key: Unique NCBI-issued identifier to enhance retrieval speed.
    :type api_key: str
    """
    if tool is not None:
        Entrez.tool = tool
    if email is not None:
        Entrez.email = email
    if api_key is not None:
        Entrez.api_key = api_key


def esearch_taxa(taxon_id, db="nucleotide", idtype="acc", retmax=99999):
    taxon_identifier_ref = "".join(["txid", str(taxon_id), "[Orgn]"])

    handle = Entrez.esearch(db=db, term=taxon_identifier_ref, idtype=idtype,
                            retmax=retmax)
    record = Entrez.read(handle)

    return record["IdList"]


def run_esearch(db="", term="", usehistory="", idtype="", retmax=99999):
    """Search for valid records in NCBI.

    Uses NCBI esearch implemented through BioPython Entrez.

    :param db: Name of the database to search.
    :type db: str
    :param term: Search term.
    :type term: str
    :param usehistory: Indicates if prior searches should be used.
    :type usehistory: str
    :return: Results of the search for each valid record.
    :rtype: dict
    """
    search_handle = Entrez.esearch(db=db, term=term, usehistory=usehistory,
                                   idtype=idtype, retmax=retmax)
    search_record = Entrez.read(search_handle)
    search_handle.close()

    return search_record


def get_records(accession_list, db="nucleotide", rettype="gb", retmode="text"):
    """Retrieve records from NCBI from a list of active accessions.

    Uses NCBI efetch implemented through BioPython Entrez.

    :param accession_list: List of NCBI accessions.
    :type accession_list: list
    :param db: Name of the database to get summaries from (e.g. 'nucleotide').
    :type db: str
    :param rettype: Type of record to retrieve (e.g. 'gb').
    :type rettype: str
    :param retmode: Format of data to retrieve (e.g. 'text').
    :type retmode: str
    :return: List of BioPython SeqRecords generated from GenBank records.
    :rtype: list
    """
    retrieved_records = []
    fetch_query = ",".join(accession_list)
    fetch_handle = Entrez.efetch(db=db, id=fetch_query, rettype=rettype,
                                 retmode=retmode)
    fetch_records = SeqIO.parse(fetch_handle, "genbank")
    for record in fetch_records:
        retrieved_records.append(record)
    fetch_handle.close()

    return retrieved_records
