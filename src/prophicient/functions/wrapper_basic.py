from subprocess import PIPE, Popen, DEVNULL

from prophicient.functions.run import run


def autoannotate(filepath, output_dir):
    """
    Annotate protein-coding genes using Prodigal

    :param filepath: path to the input FASTA file
    :type filepath: pathlib.Path
    :param output_dir: path to where the output file can be written
    :type output_dir: pathlib.Path
    """
    output_file = output_dir.joinpath(filepath.stem + ".fasta")
    command = f"prodigal -i {str(filepath)} -a {str(output_file)}"
    run(command)

    return output_file


def hhsearch(input_file, output_dir, database, cutoff):
    # take in three arguments instead
    """
    Wrapper for HHSearch
    :param args: list of parsed arguments
    :type: list
    """

    output_file = output_dir.joinpath(input_file.stem).with_suffix(".hhr")

    # run hhsearch
    # subprocess.run(["hhsearch", "-i", input_file,
    #                "-d", database, "-e", cutoff])

    # return output file instead of out

    # with Popen(["hhsearch", "-i", input_file,
    #             "-d", database, "-e", cutoff],
    #            stdout=output, stderr=output) as p:
    # out = p.stdout.read().decode("utf-8")
    # err = p.stderr.read().decode("utf-8")

    with Popen(["hhsearch", "-i", input_file, "-o", output_file,
                "-d", database, "-e", str(cutoff)],
               stdout=PIPE, stderr=DEVNULL) as process:
        process.communicate()

    return output_file