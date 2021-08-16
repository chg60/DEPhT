import shutil

from bitarray import bitarray, util as bit_util

from depht.functions.fasta import parse_fasta, write_fasta
from depht.functions.run_command import run_command

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
ENDIANESS = "big"
WORKING_DIR_NAME = "mmseqs"
FASTA_PATH_NAME = "genes"


# PHAMERATION DEFAULTS
# -----------------------------------------------------------------------------
CMODE = 0
CSTEP = 1
SENS = 8
IDENT = 0.5
COVER = 0.80
EVALUE = 0.001

MIN_GCS = 0.7


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def assemble_bacterial_mask(contigs, bacterial_fasta, gene_bit_value_path,
                            working_dir):
    """Creates a binary mask for annotated coding regions in the input
    sequence contigs by using mmseqs2 clustering to identify bacterial-looking
    genes.

    :param contigs: Input sequence contig seqrecords
    :type contigs: list[SeqRecord]
    :param bacterial_fasta: Path to a fasta of reference bacterial genes
    :type bacterial_fasta: pathlib.Path
    :param gene_bit_value_path: Path to a hexadecimal clade representation file
    :type gene_bit_value_path: pathlib.Path
    :param working_dir: Path to place created data files.
    :type working_dir: pathlib.Path
    :return: A binary bacterial mask for each inputted sequence contig
    :rtype: list[list[int]]
    """
    # Parse hexadecimal clade representation file
    bacterial_gene_bit_values, clades = parse_gene_bit_value_file(
                                                        gene_bit_value_path)

    # Create working file
    fasta_path = working_dir.joinpath(FASTA_PATH_NAME).with_suffix(".fasta")

    # Write and index input genes and initialize bacterial mask
    bacterial_masks, gene_bit_values = initialize_bacterial_mask(
                                         contigs, bacterial_fasta, fasta_path)

    # Cluster genes
    clustering_map = cluster_bacterial_genes(fasta_path, working_dir)

    # Assign a clade representation bit array to each input gene
    assign_gene_bit_values(clustering_map, bacterial_gene_bit_values,
                           gene_bit_values, bacterial_masks, clades)

    # Create a gene bit array mask from the value of the most represented clade
    clade_bit_mask = assign_clade(gene_bit_values)

    # Mark input genes well represented in the assigned clade
    mark_bacterial_mask(bacterial_masks, gene_bit_values,
                        clade_bit_mask)

    # Write bacterial masks to file
    dump_bacterial_masks(contigs, bacterial_masks, clade_bit_mask, working_dir)

    # Store bacterial mask within each contig object
    for i in range(len(contigs)):
        contig = contigs[i]
        mask = bacterial_masks[i]

        contig.update_mask_bits(mask)

    # Cleanup working file
    fasta_path.unlink()

    return bacterial_masks


def initialize_bacterial_mask(contigs, bacterial_fasta, out_path):
    """Initialize a binary bacterial mask, index, and write each gene
    in the inputted sequence contigs.

    :param contigs:  Input sequence contig seqrecords
    :type contigs: list[SeqRecord]
    :param b_gene_names: Reference bacterial gene names
    :type b_gene_names: list[str]
    :param b_gene_seqs: Reference bacterial gene sequences
    """
    # Parse reference bacterial gene multiple-sequence fasta
    b_gene_names, b_gene_seqs = parse_fasta(bacterial_fasta)
    b_gene_names = list(range(len(b_gene_names)))

    # Initialize bacterial mask data structures
    bacterial_masks = []
    gene_bit_values = []
    for contig_index, contig in enumerate(contigs):
        feature_index = 0
        for feature in contig.genes:
            # Index input gene
            b_gene_names.append(f"{contig_index}_{feature_index}")

            # Add input gene translation
            b_gene_seqs.append(str(feature.qualifiers["translation"][0]))

            feature_index += 1

        # Initialize default values for each input gene
        bacterial_masks.append([1] * feature_index)
        gene_bit_values.append([None] * feature_index)

    # Write bacterial reference genes and input genome genes to file
    write_fasta(b_gene_names, b_gene_seqs, out_path)

    return bacterial_masks, gene_bit_values


def cluster_bacterial_genes(fasta_path, working_dir):
    """Clusters input genes with MMseqs2.

    :param fasta_path: Path to a fasta with input and reference genes.
    :type fasta_path: pathlib.Path
    :param working_dir: Path to place created data files
    :type working_dir: pathlib.Path
    """
    # Create database file directory
    database_dir = working_dir.joinpath("database")
    database_dir.mkdir(exist_ok=True)

    # Define mmseqs working file paths
    sequence_db = database_dir.joinpath("sequenceDB")
    cluster_db = database_dir.joinpath("clusterDB")
    seqfile_db = database_dir.joinpath("seqfileDB")
    result_file = database_dir.joinpath("clustering_results.txt")
    tmp_dir = database_dir.joinpath("tmp")
    tmp_dir.mkdir(exist_ok=True)

    # Create the initial mmseqs DB
    mmseqs_createdb(fasta_path, sequence_db)

    # Cluster bacterial genes in linear time
    mmseqs_linclust(sequence_db, cluster_db, tmp_dir, CMODE,
                    IDENT, COVER, EVALUE)

    # Create a database from the clustering
    mmseqs_createseqfiledb(sequence_db, cluster_db, seqfile_db)

    # Write the clustering results to file
    mmseqs_result2flat(sequence_db, sequence_db, seqfile_db,
                       result_file)

    # Parse the clustering output
    clustering_results = parse_mmseqs(result_file)

    # Remove database file directory
    shutil.rmtree(database_dir)

    return clustering_results


def assign_gene_bit_values(clustering_map, bacterial_gene_bit_values,
                           gene_bit_values, bacterial_masks, clades):
    """Assigns genes values contained in bit arrays depending on clade
    representation of reference bacterial genes.

    :param clustering_map: Cluster ID to cluster gene-member map
    :type clustering_map: dict
    :param bacterial_gene_bit_values: Bit arrays for bacterial reference genes.
    :type bacterial_gene_bit_values: list
    :param gene_bit_values: Initialized gene bit value list
    :type gene_bit_values: list
    :param bacterial_masks: Initialized bacterial core gene bits
    :type bacterial_masks: list
    :param clades: Number of total clades
    :type clades: int
    """
    for cluster, gene_names in clustering_map.items():
        gois = []
        cluster_bitarray = None
        absolute_bitarray = None
        for gene_name in gene_names:
            # Identify and store indexed input sequence genes
            if "_" in gene_name:
                gois.append(gene_name)
                continue

            # Reidentify bacterial reference gene index
            gene_index = int(gene_name)

            # Initialize / OR cluster member bitarray(s)
            # given reference gene values, to determine clade representation
            if cluster_bitarray is None:
                cluster_bitarray = bacterial_gene_bit_values[gene_index]
            else:
                # Ensure input bitarrays are the same length
                cluster_bitarray, b_gene_bitarray = equalize_bitarrays(
                                        cluster_bitarray,
                                        bacterial_gene_bit_values[gene_index])

                # OR reference gene bitarrays
                cluster_bitarray = (cluster_bitarray | b_gene_bitarray)

            # Intialize / AND cluster member bitarrays()
            # given reference gene values, to determine high confidence
            # bacterial genes
            if absolute_bitarray is None:
                absolute_bitarray = bacterial_gene_bit_values[gene_index]
            else:
                absolute_bitarray, b_gene_bitarray = equalize_bitarrays(
                                        absolute_bitarray,
                                        bacterial_gene_bit_values[gene_index])

                absolute_bitarray = (absolute_bitarray & b_gene_bitarray)

        for goi in gois:
            goi_split = goi.split("_")

            # Retrieve contig and gene indicies from input gene name
            contig_i = int(goi_split[0])
            gene_i = int(goi_split[1])

            # Assign input gene a value based on it's reference relatives
            gene_bit_values[contig_i][gene_i] = cluster_bitarray

            if absolute_bitarray is not None:
                if absolute_bitarray.count() >= clades:
                    bacterial_masks[contig_i][gene_i] = 0


def assign_clade(gene_bit_values, min_gcs=MIN_GCS):
    """Assigns clade membership of the input sequences based on the majority
    gene bit value

    :param gene_bit_values: Bitarrays of the clade representation for each gene
    :type gene_bit_values: list[bitarray]
    :return: A bitarray mask for the majority-assigned clade
    :rtype: bitarray
    """
    bit_count = list()
    for contig_gene_bitarrays in gene_bit_values:
        for gene_bitarray in contig_gene_bitarrays:
            if gene_bitarray is None:
                continue

            for index in range(len(gene_bitarray)):
                bit = gene_bitarray[index]
                if index >= len(bit_count):
                    bit_count.append(0)

                if bit:
                    bit_count[index] += 1

    if len(bit_count) <= 0:
        return bitarray([0])

    count = max(bit_count)
    if (count / len(gene_bit_values)) < min_gcs:
        return bitarray([0] * len(bit_count), endian=ENDIANESS)

    clade = bit_count.index(max(bit_count))

    clade_bit_mask = bitarray([0] * len(bit_count), endian=ENDIANESS)
    clade_bit_mask[clade] = 1

    return clade_bit_mask


def mark_bacterial_mask(bacterial_masks, gene_bit_values, clade_bit_mask):
    """Marks genes as bacterial with the given clade value bit mask.

    :param bacterial_masks: Initialized binary bacterial mask list
    :type bacterial_masks: list[list[int]]
    :param gene_bit_values: Bitarrays of the clade representation for each gene
    :type gene_bit_values: list[bitarray]
    :param clade_bit_mask: Bitarray mask for the majority-asssigned clade
    :type clade_bit_mask: bitarray
    """
    for contig_i, bacterial_mask in enumerate(bacterial_masks):
        for gene_i, bacterial_bit in enumerate(bacterial_mask):
            if bacterial_bit == 0:
                continue

            gene_bitarray = gene_bit_values[contig_i][gene_i]

            if gene_bitarray is None:
                continue

            masked_bitarray = (gene_bitarray & clade_bit_mask)

            if masked_bitarray.count() > 0:
                bacterial_mask[gene_i] = 0


def dump_bacterial_masks(contigs, bacterial_masks, clade_bit_mask,
                         working_dir):
    for i, contig in enumerate(contigs):
        bacterial_mask = bacterial_masks[i]
        filepath = working_dir.joinpath(".".join([contig.id, "txt"]))

        filehandle = filepath.open(mode="w")

        for bit in clade_bit_mask:
            filehandle.write(str(int(bit)))

        filehandle.write("\n")

        for bit in bacterial_mask:
            filehandle.write(str(int(bit)))


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def equalize_bitarrays(a_bitarray, b_bitarray):
    length_disparity = len(a_bitarray) - len(b_bitarray)

    if length_disparity < 0:
        for _ in range(abs(length_disparity)):
            a_bitarray.append(0)
    elif length_disparity > 0:
        for _ in range(length_disparity):
            b_bitarray.append(0)

    return a_bitarray, b_bitarray


# MMSEQS HELPER FUNCTIONS
# -----------------------------------------------------------------------------

def mmseqs_createdb(fasta, mmseqsdb):
    """
    Creates an MMseqs2 database at `mmseqsdb` using `fasta` as its input.

    :param fasta: path to the FASTA file to assemble into an MMseqs2 database
    :type fasta: pathlib.Path
    :param mmseqsdb: path to the desired MMseqs2 database
    :type mmseqsdb: pathlib.Path
    :return:
    """
    c = f"mmseqs createdb {fasta} {mmseqsdb} -v 3"
    run_command(c)


def mmseqs_cluster(sequence_db, cluster_db, tmp_dir, clustermode, clustersteps,
                   sensitivity, minseqid, coverage, evalue):
    """
    Clusters an MMseqs2 sequence database using the given parameters.

    :param sequence_db: path to the MMseqs2 sequence database to cluster
    :type sequence_db: pathlib.Path
    :param cluster_db: path to the MMseqs2 output cluster database
    :type cluster_db: pathlib.Path
    :param tmp_dir: temporary directory that MMseqs2 can use
    :type tmp_dir: pathlib.Path
    :param clustermode: mmseqs cluster --cluster-mode
    :type clustermode: int
    :param clustersteps: mmseqs cluster --cluster-steps
    :type clustersteps: int
    :param sensitivity: mmseqs cluster -s
    :type sensitivity: float
    :param minseqid: mmseqs cluster --min-seq-id
    :type minseqid: float
    :param coverage: mmseqs cluster -c
    :type coverage: float
    :param evalue: mmseqs cluster -e
    :type evalue: float
    :return:
    """
    c = (f"mmseqs cluster {str(sequence_db)} {str(cluster_db)} {str(tmp_dir)} "
         f"-v 3 --max-seqs 1000 --cluster-mode {clustermode} "
         f"--cluster-steps {clustersteps} -s {sensitivity} "
         f"--min-seq-id {minseqid} -c {coverage} -e {evalue}")
    run_command(c)


def mmseqs_linclust(sequence_db, cluster_db, tmp_dir, clustermode,
                    minseqid, coverage, evalue):
    """
    Clusters an MMseqs2 sequence database using the given parameters in
    linear time.

    :param sequence_db: path to the MMseqs2 sequence database to cluster
    :type sequence_db: pathlib.Path
    :param cluster_db: path to the MMseqs2 output cluster database
    :type cluster_db: pathlib.Path
    :param tmp_dir: temporary directory that MMseqs2 can use
    :type tmp_dir: pathlib.Path
    :param clustermode: mmseqs cluster --cluster-mode
    :type clustermode: int
    :param minseqid: mmseqs cluster --min-seq-id
    :type minseqid: float
    :param coverage: mmseqs cluster -c
    :type coverage: float
    :param evalue: mmseqs cluster -e
    :type evalue: float
    :return:
    """
    c = (f"mmseqs linclust {sequence_db} {cluster_db} {tmp_dir} "
         f"-v 3 --cluster-mode {clustermode} "
         f"--min-seq-id {minseqid} -c {coverage} -e {evalue}")
    run_command(c)


def mmseqs_result2profile(sequence_db, cluster_db, profile_db):
    """
    Converts an MMseqs2 cluster output database to a profile database.

    :param sequence_db: path to the clustered MMseqs2 sequence database
    :type sequence_db: pathlib.Path
    :param cluster_db: path to the MMseqs2 cluster database
    :type cluster_db: pathlib.Path
    :param profile_db: path to the desired MMseqs2 profile database
    :type profile_db: pathlib.Path
    :return:
    """
    c = (f"mmseqs result2profile {str(sequence_db)} {str(sequence_db)} "
         f"{str(cluster_db)} {str(profile_db)} -v 3")
    run_command(c)


def mmseqs_profile2consensus(profile_db, consensus_db):
    """
    Extracts consensus sequences from an MMseqs2 profile database and creates
    an MMseqs2 sequence database from those consensus sequences.

    :param profile_db: path to the MMseqs2 profile database
    :type profile_db: pathlib.Path
    :param consensus_db: path to a MMseqs2 consensus sequence database
    :type consensus_db: pathlib.Path
    :return:
    """
    c = f"mmseqs profile2consensus {profile_db} {consensus_db} -v 3"
    run_command(c)


def mmseqs_search(profile_db, consensus_db, align_db, tmp_dir, minseqid,
                  coverage, evalue):
    """
    Searches an MMseqs2 profile database against an MMseqs2 profile consensus
    sequence database, for HMM-based clustering.

    :param profile_db: path to the MMseqs2 profile database
    :type profile_db: pathlib.Path
    :param consensus_db: path to the MMseqs2 consensus sequence database
    :type consensus_db: pathlib.Path
    :param align_db: path to the desired MMseqs2 search result database
    :type align_db: pathlib.Path
    :param tmp_dir: temporary directory that MMseqs2 can use
    :type tmp_dir: pathlib.Path
    :param minseqid: mmseqs search --min-seq-id
    :type minseqid: float
    :param coverage: mmseqs search -c
    :type coverage: float
    :param evalue: mmseqs search --e-profile
    :type evalue: float
    :return:
    """
    c = (f"mmseqs search {profile_db} {consensus_db} {align_db} {tmp_dir} "
         f"-v 3 --max-seqs 1000 --min-seq-id {minseqid} -c {coverage} "
         f"--cov {coverage} -e {evalue} --e-profile {evalue} "
         "--add-self-matches")
    run_command(c)


def mmseqs_clust(consensus_db, align_db, result_db):
    """
    Clusters the MMseqs2 result database from running 'mmseqs search'.

    :param consensus_db: path to the MMseqs2 consensus sequence database
    :type consensus_db: pathlib.Path
    :param align_db: path to the MMseqs2 search result database
    :type align_db: pathlib.Path
    :param result_db: path to the desired MMseqs2 clustered result database
    :type result_db: pathlib.Path
    :return:
    """
    c = f"mmseqs clust {consensus_db} {align_db} {result_db} -v 3"
    run_command(c)


def mmseqs_createseqfiledb(sequence_db, cluster_db, sf_db):
    """
    Creates an MMseqs2 seqfileDB from a clustering output database.

    :param sequence_db: the path to the query db in MMseqs2 clustering
    :type sequence_db: pathlib.Path
    :param cluster_db: the path to the MMseqs2 clustering result database
    :type cluster_db: pathlib.Path
    :param sf_db: the path to the desired seqfileDB
    :type sf_db: pathlib.Path
    :return:
    """
    c = f"mmseqs createseqfiledb {sequence_db} {cluster_db} {sf_db} -v 3"
    run_command(c)


def mmseqs_result2flat(query_db, subject_db, result_db, output):
    """
    Converts an MMseqs2 result seqfileDB to a FASTA-like output file.

    :param query_db: the database used as the query in MMseqs2 clustering
    :type query_db: pathlib.Path
    :param subject_db: the database used as the subject in MMseqs2 clustering
    :type subject_db: pathlib.Path
    :param result_db: the MMseqs2 clustering result database
    :type result_db: pathlib.Path
    :param output: the desired clustering output file
    :type output: pathlib.Path
    :return:
    """
    c = f"mmseqs result2flat {query_db} {subject_db} {result_db} {output} -v 3"
    run_command(c)


# FILEIO
# -----------------------------------------------------------------------------
def parse_mmseqs(filepath):
    """
    Parses the indicated MMseqs2 output into a dictionary mapping
    phamids to the geneids that were put in them

    :param filepath: path to the MMseqs2 output file to parse
    :type filepath: pathlib.Path
    :return: phams
    """
    phams = dict()
    phamid = 0
    pham_genes = list()

    with filepath.open("r") as fh:
        prior = fh.readline()
        current = fh.readline()

        # While loop to iterate until EOF
        while current:
            if current.startswith(">"):
                # If current & prior both header lines - new pham block
                if prior.startswith(">"):
                    try:
                        pham_genes.pop(-1)
                    except IndexError:
                        pass
                    phams[phamid] = pham_genes
                    phamid += 1
                    pham_genes = list()
                pham_genes.append(current.lstrip(">").rstrip())
            else:
                # Do nothing on translation lines
                pass
            prior, current = current, fh.readline()
        # Need to dump the last p into the dictionary
        phams[phamid] = pham_genes
    phams.pop(0)  # 0th pham is placeholder

    return phams


def parse_gene_bit_value_file(filepath):
    """
    Parses the gene bit value file for determining input genome clade
    membership

    :param filepath: path to the gene bit value file
    :type filepath: pathlib.Path
    :return: Returns a list of bit arrays
    :rtype: list[bitarray]
    """
    bacterial_gene_bit_values = []
    filehandle = open(filepath, "rb")
    line = filehandle.readline()

    hex_values = line.split(b"_")

    clades = 0
    for hex_value in hex_values:
        if hex_value == b"":
            continue

        b_gene_bitarray = bit_util.hex2ba(hex_value, endian=ENDIANESS)

        bitarray_count = b_gene_bitarray.count()
        if bitarray_count > clades:
            clades = bitarray_count

        bacterial_gene_bit_values.append(b_gene_bitarray)

    filehandle.close()

    return bacterial_gene_bit_values, clades
