from depht.functions.subprocess import run_command


def clustalo(fasta_file, aln_out_path, mat_out_path=None, tree_out_path=None,
             outfmt="fasta", infmt="fasta", threads=1, verbose=0):
    """
    Runs Clustal Omega to generate a multiple sequence alignment (MSA)
    and percent identity matrix (PIM) for the indicated file. Infile is
    expected to be in FASTA multiple sequence format. MSA will be in
    Clustal format.
    :param fasta_file: FASTA file containing sequences to be aligned
    :type fasta_file: str or pathlib.Path
    :param aln_out_path: the multiple sequence alignment (MSA) output file
    :type aln_out_path: str or pathlib.Path
    :param mat_out_path: The percent identity matrix (PIM) output file
    :type mat_out_path: str or pathlib.Path
    :param tree_out_path: The alignment guide tree output file
    :type tree_out_path: str or pathlib.Path
    :param outfmt: The file format of the alignment to be exported.
    :type outfmt: str
    :param infmt:  The file format of the sequence file to be read in
    :type infmt: str
    :param threads: number of threads to use for alignment
    :type threads: int
    :param verbose: verbosity level (0-2)
    :type verbose: int
    :return: Returns the sequence alignment path[0] and distance matrix path[1]
    :rtype: tuple
    """
    # Make sure verbose is in range(0,3)
    if verbose <= 0:
        verbose = 0
    elif verbose > 2:
        verbose = 2

    # Build Clustal Omega command that will produce a clustal-formatted
    # alignment output file and percent identity matrix
    command = f"""clustalo -i "{fasta_file}" -o "{aln_out_path}" """ \
              f" --outfmt={outfmt} --infmt={infmt} "\
              f"--force --output-order=tree-order " \
              f"--threads={threads}"

    if mat_out_path is not None:
        command = " ".join([command, (f"""--distmat-out="{mat_out_path}" """
                                      "--full --percent-id")])
    if tree_out_path is not None:
        command = " ".join([command, (
                                f"""--guidetree-out="{tree_out_path}" """)])

    for _ in range(verbose):
        command += " -v"                # Add verbosity to command

    # Run the Clustal Omega command as a subprocess
    run_command(command)

    # Return output files so user can find them easily, whether they specified
    # the filenames or not
    return (aln_out_path, mat_out_path)
