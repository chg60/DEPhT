"""Utility script to perform pham assembly with MMseqs2.
"""
import argparse
import json
import pathlib
import shutil
import sys

from depht.functions import mmseqs
from depht.functions.fasta import parse_fasta
from depht_train.classes.database import Database, Pham


def parse_args(unparsed_args):
    """Parse commandline arguments.

    :param unparsed_args:
    :type unparsed_args: list[str]
    """
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("fasta", type=pathlib.Path,
                   help="path to the FASTA file containing genes to phamerate")
    p.add_argument("outdir", type=pathlib.Path,
                   help="directory where file I/O can occur")
    p.add_argument("param_file", type=pathlib.Path,
                   help="path to the JSON file of MMseqs2 parameters to use")

    return p.parse_args(unparsed_args)


def parse_param_file(param_file):
    """Attempts to parse a given parameter JSON file to two
    tuples of first-iter and second-iter parameter sets for MMseqs2.

    :param param_file: the JSON-formatted parameter file
    :type param_file: pathlib.Path
    :return: first_iter, second_iter
    """
    with param_file.open("r") as param_handle:
        param_dict = json.load(param_handle)

    return parse_param_dict(param_dict)


def parse_param_dict(param_dict):
    """
    Attempts to parse a given parameter dictionary to two
    tuples of first-iter and second-iter parameter sets for MMseqs2.

    :param param_dict: the formatted parameter dictionary
    :type param_dict: dict
    :return: first_iter, second_iter
    """
    first_iter_dict = param_dict.get("first_iteration")
    a = first_iter_dict.get("--cluster-mode")       # Define clustermode
    b = first_iter_dict.get("--cluster-steps")      # Define clustersteps
    c = first_iter_dict.get("-s")                   # Define sensitivity
    d = first_iter_dict.get("--min-seq-id")         # Define min-seq-id
    e = first_iter_dict.get("-c")                   # Define coverage
    f = first_iter_dict.get("-e")                   # Define e-value
    first_iter = (a, b, c, d, e, f)

    second_iter_dict = param_dict.get("second_iteration")
    if second_iter_dict is None:
        second_iter = None
    else:
        g = second_iter_dict.get("--min-seq-id")        # Define min-seq-id
        h = second_iter_dict.get("-c")                  # Define coverage
        j = second_iter_dict.get("-e")                  # Define e-value
        second_iter = (g, h, j)

    return first_iter, second_iter


def phamerate(sequence_db, db, tmpdir, first_iter, second_iter=None):
    """

    :param sequence_db: path to the MMseqs2 sequence database
    :type sequence_db: pathlib.Path
    :param db: database of geneids/translations being phamerated
    :type db: Database
    :param tmpdir:
    :type tmpdir: pathlib.Path
    :param first_iter:
    :type first_iter: tuple
    :param second_iter:
    :type second_iter: tuple
    :return:
    """
    # Refresh tmpdir
    if not tmpdir.is_dir():
        tmpdir.mkdir(parents=False)

    temp_phams = list()

    # First iteration paths
    cluster_db = tmpdir.joinpath("clusterDB")
    seqfile_db = tmpdir.joinpath("seqfileDB")
    first_iter_out = tmpdir.joinpath("first_iter_out.txt")

    # Expand first iteration parameter set
    cmode, cstep, sens, ident, cover, evalue = first_iter

    # Cluster database using given parameters
    mmseqs.mmseqs_cluster(sequence_db, cluster_db, tmpdir, cmode, cstep, sens,
                          ident, cover, evalue)
    mmseqs.mmseqs_createseqfiledb(sequence_db, cluster_db, seqfile_db)
    mmseqs.mmseqs_result2flat(sequence_db, sequence_db, seqfile_db,
                              first_iter_out)
    first_iter_phams = mmseqs.parse_mmseqs(first_iter_out)

    if second_iter:
        # Second iteration paths
        profile_db = tmpdir.joinpath("profileDB")
        consensus_db = tmpdir.joinpath("consensusDB")
        align_db = tmpdir.joinpath("alignDB")
        result_db = tmpdir.joinpath("resultDB")
        hmm_seqfile_db = tmpdir.joinpath("hmmSeqfileDB")
        second_iter_out = tmpdir.joinpath("second_iter_out.txt")

        # Expand second iteration parameter set
        ident, cover, evalue = second_iter

        # Cluster profiles using given parameters
        mmseqs.mmseqs_result2profile(sequence_db, cluster_db, profile_db)
        mmseqs.mmseqs_profile2consensus(profile_db, consensus_db)
        mmseqs.mmseqs_search(profile_db, consensus_db, align_db, tmpdir,
                             ident, cover, evalue)
        mmseqs.mmseqs_clust(consensus_db, align_db, result_db)
        mmseqs.mmseqs_createseqfiledb(sequence_db, result_db, hmm_seqfile_db)
        mmseqs.mmseqs_result2flat(sequence_db, sequence_db, hmm_seqfile_db,
                                  second_iter_out)
        second_iter_phams = mmseqs.parse_mmseqs(second_iter_out)

        lookup = dict()
        for phamid, pham_geneids in first_iter_phams.items():
            for pham_geneid in pham_geneids:
                lookup[pham_geneid] = phamid

        for source_id, source_geneids in second_iter_phams.items():
            all_pham_geneids = set()
            for source_geneid in source_geneids:
                target_id = lookup[source_geneid]
                for target_geneid in first_iter_phams[target_id]:
                    all_pham_geneids.add(target_geneid)
            all_pham_geneids = list(all_pham_geneids)
            all_pham_translations = [db.get_translation_from_geneid(x) for x in
                                     all_pham_geneids]
            temp_phams.append(Pham(all_pham_geneids, all_pham_translations))

    else:
        for pham_id, all_pham_geneids in first_iter_phams.items():
            all_pham_translations = [db.get_translation_from_geneid(x) for x in
                                     all_pham_geneids]
            temp_phams.append(Pham(all_pham_geneids, all_pham_translations))

    return temp_phams


def execute_phamerate_pipeline(fasta, outdir, first_iter, second_iter=None):
    tempdir = outdir.joinpath("temp")

    # Make sure outdir and tmpdir exist
    if not outdir.is_dir():
        outdir.mkdir(parents=True)
    if not tempdir.is_dir():
        tempdir.mkdir(parents=False)

    # Read fasta file and build database instance
    geneids, translations = parse_fasta(fasta)
    input_db = Database(geneids, translations)

    # Create mmseqsdb
    nr_fasta = outdir.joinpath("genes.fasta")
    with nr_fasta.open("w") as fh:
        fh.write(repr(input_db))
    mmseqsdb = tempdir.joinpath("sequenceDB")
    mmseqs.mmseqs_createdb(nr_fasta, mmseqsdb)

    # Phamerate!
    phams = phamerate(mmseqsdb, input_db, tempdir, first_iter,
                      second_iter)

    # Write fasta file for each pham
    for i, pham in enumerate(phams):
        pham_fasta = outdir.joinpath(f"pham_{i}.fasta")
        pham_translations = pham.get_translations()

        for translation in pham_translations:
            for geneid in input_db.get_geneids_from_translation(translation):
                try:
                    pham.add_gene(geneid, translation)
                except ValueError:
                    pass    # geneid already in

        with pham_fasta.open("w") as fh:
            fh.write(str(pham))

    # Clean up temporary directory
    shutil.rmtree(tempdir)

    # Clean up temporary mulitple sequence fasta file
    if nr_fasta.is_file():
        nr_fasta.unlink()


def main(unparsed_args=None):
    """Commandline entrypoint to this module."""
    if not unparsed_args:
        unparsed_args = sys.argv

    if len(unparsed_args) == 1:
        unparsed_args.append("-h")

    args = parse_args(unparsed_args[1:])

    first_iter_params, second_iter_params = parse_param_file(args.param_file)

    execute_phamerate_pipeline(args.fasta, args.outdir,
                               first_iter_params, second_iter_params)


if __name__ == "__main__":
    main()
