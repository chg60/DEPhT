import subprocess


def prodigal(args):
    """
    Wrapper function for prodigal
    :param args: list of arguments
    :type args: list
    """
    # the formatted arguments to run prodigal will be stored in a list
    formatted_args = ["prodigal"]

    # names the output file phage_name.output_format
    phage_name = args.input_path.stem

    output_format = ".gbk"  # default

    output_filename = phage_name + output_format

    # creating the output filepath
    output_filepath = args.output_path/output_filename

    formatted_args.extend(["-i", args.input_path, "-a", output_filepath])

    if args.training is not None:
        formatted_args.extend(["-t", args.training])
    if args.nucleotide is not None:
        formatted_args.extend(["-d", args.nucleotide])
    if args.genes is not None:
        formatted_args.extend(["-s", args.genes])

    subprocess.run(formatted_args)

    return output_filepath


def hhsearch(input_file, output_dir, database, cutoff):
   # take in three arguments instead
    """
    Wrapper for HHSearch
    :param args: list of parsed arguments
    :type: list
    """

    output_file = Path(output_dir/(input_file.stem + ".hhr"))

    print(output_file)

    with open(output_file, "w") as output:

        # run hhsearch
        # subprocess.run(["hhsearch", "-i", input_file, "-d", database, "-e", cutoff])

        # return output file instead of out

        # with Popen(["hhsearch", "-i", input_file, "-d", database, "-e", cutoff], stdout=output, stderr=output) as p:
        # out = p.stdout.read().decode("utf-8")
        #err = p.stderr.read().decode("utf-8")

        with Popen(["hhsearch", "-i", input_file, "-d", database, "-e", cutoff],
                stdout=output, stderr=output) as process:
            process.communicate()

    return output_file
