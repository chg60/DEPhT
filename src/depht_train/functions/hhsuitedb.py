import json
import os
import shlex
import shutil
from subprocess import DEVNULL, PIPE, Popen


# CREATE DATABASE
# -----------------------------------------------------------------------------
def create_hhsuitedb(aln_dir, db_dir, db_name, cores=1, use_mpi=False,
                     verbose=False, versions=None):
    """Function to create a database compatible with the HHsuite3 toolkit.

    :param aln_dir: Path to a directory containing fasta-formatted MSA files.
    :type aln_dir: pathlib.Path
    :param db_dir: Path to a directory where the database will be created.
    :type db_dir: pathlib.Path
    :param db_name: Name of the database to be created.
    :type db_name: str
    :param cores: Number of subprocesses to spawn during database construction.
    :type cores: int
    :param use_mpi: Toggles use of the mpirun library for parallel-processing.
    :type use_mpi: bool
    """
    if verbose:
        stdout = PIPE
    else:
        stdout = DEVNULL

    if verbose:
        print("Creating multiple sequence ffindex file...")
    msa_fftuple = build_msa_ffindex(aln_dir, db_dir, db_name, stdout=stdout)

    # A hhsuite database requires pairs of three main file types:
    # A .ffdata and .ffindex pair containing a3m formatted MSAs
    if verbose:
        print("Condensing multiple sequence ffindex file into a3m file...")
    a3m_fftuple = convert_a3m_ffindex(db_dir, db_name, msa_fftuple,
                                      cores=cores, stdout=stdout,
                                      use_mpi=use_mpi)

    # A .ffdata and .ffindex pair containing formatted HHsuite3 HMMs
    if verbose:
        print("Creating Hidden Markov Model ffindex file...")
    hmm_fftuple = convert_hmm_ffindex(db_dir, db_name, a3m_fftuple,
                                      cores=cores, stdout=stdout,
                                      use_mpi=use_mpi)

    # A .ffdata and .ffindex pair containing cs219 sequences for prefiltering
    if verbose:
        print("Creating cs219 ffindex file...")
    cs219_ffindex = create_cs219_ffindex(db_dir, db_name, stdout=stdout)

    # The MSA and HMM ffindex files are sorted according to data from cs219
    # These steps are essential to the speed of the HHsuite3 utilities
    if verbose:
        print("Sorting ffindex files for fast database indexing...")
    sorting_file = create_sorting_file(cs219_ffindex, db_dir, stdout=stdout)
    sort_ffindex(db_dir, sorting_file, a3m_fftuple, stdout=stdout)
    sort_ffindex(db_dir, sorting_file, hmm_fftuple, stdout=stdout)

    # Removes msa_ffindex files, as they are not needed for a functional db
    msa_fftuple[0].unlink()
    msa_fftuple[1].unlink()
    sorting_file.unlink()

    # Is a placeholder.
    # Requires additional work for the functionality to be realized.
    if not verify_hhsuite_database(db_dir, db_name):
        print(f"Inconsistencies detected in HHsuite database "
              f"at '{db_dir}'.\n  Scrapping database.")
        shutil.rmtree(db_dir)

    # Creates a json file that dumps key: value pairs with no file extension
    # Is unnecessesary, only for convenience and easy tab completion
    db_path = db_dir.joinpath(db_name)
    if versions is not None:
        create_version_file(db_path, versions)

    return db_path


def build_msa_ffindex(aln_dir, db_dir, db_name, stdout=None):
    """Function to build a ffindex database pair from a directory of
    fasta-formatted MSA files

    :param aln_dir: Path to a directory containing fasta-formatted MSA files.
    :type aln_dir: pathlib.Path
    :param db_dir: Path to a directory where the database will be created.
    :type db_dir: pathlib.Path
    :param db_name: Name of the database to be created.
    :type db_name: str
    :param stdout: Standard output to use for the spawned process
    :returns: A tuple containing path objects for the created ffindex database
    :rtype: tuple(pathlib.Path, pathlib.Path)
    """
    if stdout is None:
        stdout = DEVNULL

    msa_ffdata = db_dir.joinpath("".join([db_name, "_msa.ffdata"]))
    msa_ffindex = db_dir.joinpath("".join([db_name, "_msa.ffindex"]))

    command = f"ffindex_build -s '{msa_ffdata}' '{msa_ffindex}' '{aln_dir}'"

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (msa_ffdata, msa_ffindex)


def build_hmm_ffindex(hmm_dir, db_dir, db_name, stdout=None):
    """Function to build a ffindex database pair from a directory of
    HHsuite3 formatted HMM files

    :param hmm_dir: Path to a directory containing HHsuite3 HMM files.
    :type aln_dir: pathlib.Path
    :param db_dir: Path to a directory where the database will be created.
    :type db_dir: pathlib.Path
    :param db_name: Name of the database to be created.
    :type db_name: str
    :param stdout: Standard output to use for the spawned process
    :returns: A tuple containing path objects for the created ffindex database
    :rtype: tuple(pathlib.Path, pathlib.Path)
    """
    if stdout is None:
        stdout = DEVNULL

    hmm_ffdata = db_dir.joinpath("".join([db_name, "_hmm.ffdata"]))
    hmm_ffindex = db_dir.joinpath("".join([db_name, "_hmm.ffindex"]))

    command = f"ffindex_build -s '{hmm_ffdata}' '{hmm_ffindex}' '{hmm_dir}'"

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (hmm_ffdata, hmm_ffindex)


def convert_a3m_ffindex(db_dir, db_name, msa_fftuple, cores=1, stdout=None,
                        use_mpi=False):
    """Function to convert a fasta-formatted MSA ffindex database into an
    a3m formatted ffindex database

    :param db_dir: Path to a directory where the database will be created.
    :type db_dir: pathlib.Path
    :param db_name: Name of the database to be created.
    :type db_name: str
    :param msa_fftuple: A tuple containing path objects to an MSA ffindex db.
    :type msa_fftuple: tuple(pathlib.Path, pathlib.Path)
    :param cores: Number of subprocesses to spawn during database construction.
    :type cores: int
    :param stdout: Standard output to use for the spawned process(es)
    :param use_mpi: Toggles use of the mpirun library for parallel-processing.
    :type use_mpi: bool
    :returns: A tuple containing path objects for the created ffindex database.
    :rtype: tuple(pathlib.Path, pathlib.Path)
    """
    if stdout is None:
        stdout = DEVNULL

    a3m_ffdata = db_dir.joinpath("".join([db_name, "_a3m.ffdata"]))
    a3m_ffindex = db_dir.joinpath("".join([db_name, "_a3m.ffindex"]))

    # This command assumes that the msa_fftuple holds the path to the ffindex
    # data file at index 0 and the ffindex index file at index 1.
    # The parameter set is supplied by the HHsuite3 documentation.
    command = (f"'{msa_fftuple[0]}' '{msa_fftuple[1]}' "
               f"-d '{a3m_ffdata}' -i '{a3m_ffindex}' "
               "-- hhconsensus -M 50 -maxres 65535 "
               "-i stdin -oa3m stdout -v 0")

    # If parallel processing with mpirun is desired, use ffindex_apply_mpi.
    # ffindex_apply_mpi is tricky to install due to dependancies, be careful.
    if use_mpi:
        os.environ["OMP_NUM_THREADS"] = "1"
        command = "".join([(f"mpirun -np {cores} ffindex_apply_mpi "),
                           command])
    else:
        command = "".join(["ffindex_apply ", command])

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (a3m_ffdata, a3m_ffindex)


def convert_hmm_ffindex(db_dir, db_name, a3m_fftuple, cores=1, stdout=None,
                        use_mpi=False):
    """Function to convert a fasta-formatted MSA ffindex database into an
    HHsuite3 formatted HMM ffindex database

    :param db_dir: Path to a directory where the database will be created.
    :type db_dir: pathlib.Path
    :param db_name: Name of the database to be created.
    :type db_name: str
    :param a3m_fftuple: A tuple containing path objects to an MSA ffindex db.
    :type a3m_fftuple: tuple(pathlib.Path, pathlib.Path)
    :param cores: Number of subprocesses to spawn during database construction.
    :type cores: int
    :param stdout: Standard output to use for the spawned process(es)
    :param use_mpi: Toggles use of the mpirun library for parallel-processing.
    :type use_mpi: bool
    :returns: A tuple containing path objects for the created ffindex database.
    :rtype: tuple(pathlib.Path, pathlib.Path)
    """
    if stdout is None:
        stdout = DEVNULL

    hmm_ffdata = db_dir.joinpath("".join([db_name, "_hmm.ffdata"]))
    hmm_ffindex = db_dir.joinpath("".join([db_name, "_hmm.ffindex"]))

    # This command assumes that the a3m_fftuple holds the path to the ffindex
    # data file at index 0 and the ffindex index file at index 1.
    command = (f"'{a3m_fftuple[0]}' '{a3m_fftuple[1]}' "
               f"-i '{hmm_ffindex}' -d '{hmm_ffdata}' "
               "-- hhmake -i stdin -o stdout -v 0")

    # If parallel processing with mpirun is desired, use ffindex_apply_mpi.
    # ffindex_apply_mpi is tricky to install due to dependancies, be careful.
    if use_mpi:
        os.environ["OMP_NUM_THREADS"] = "1"
        command = "".join([f"mpirun -np {cores} ffindex_apply_mpi ",
                           command])
    else:
        command = "".join(["ffindex_apply ", command])

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (hmm_ffdata, hmm_ffindex)


def create_cs219_ffindex(db_dir, db_name, cores=1, use_mpi=False,
                         stdout=None):
    """Function to convert a a3m formatted MSA ffindex database into
    a cs219 sequence database for fast prefiltering

    :param db_dir: Path to a directory where the database will be created.
    :type db_dir: pathlib.Path
    :param db_name: Name of the database to be created.
    :type db_name: str
    :param cores: Number of subprocesses to spawn during database construction.
    :type cores: int
    :param use_mpi: Toggles use of the mpirun library to parallel-processing.
    :type use_mpi: bool
    :param stdout: Standarad output to use for the spawned process(es)
    :returns: Path to the created cs219 sequence ffindex index file
    :rtype: pathlib.Path
    """
    if stdout is None:
        stdout = DEVNULL

    # Annoyingly, this function requires that the path given has no file ext.
    a3m_ffpath = db_dir.joinpath("".join([db_name, "_a3m"]))
    cs219_ffpath = db_dir.joinpath("".join([db_name, "_cs219"]))

    # The parameter set is supplied by the HHsuite3 documentation
    command = ("cstranslate -f -x 0.3 -c 4 -I a3m "
               f"-i '{a3m_ffpath}' -o '{cs219_ffpath}'")

    if use_mpi:
        command = "".join([f"mpirun -np {cores}", command])

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    cs219_ffindex = db_dir.joinpath("".join([db_name, "_cs219.ffindex"]))
    return cs219_ffindex


def create_sorting_file(cs219_ffindex, db_dir, sorting_file=None,
                        stdout=None):
    """Function to create a sorted index from a cs219 sequences file.

    :param cs219_ffindex: Path to a ffindex index file for cs219 sequences.
    :type cs219_ffindex: pathlib.Path
    :param db_dir: Path to a directory where the database will be created.
    :type db_dir: pathlib.Path
    :param sorting_file: Desired path to create the sorting file
    :type sorting_file: pathlib.Path
    :param stdout: Standard output to use for the spawned process
    :returns: Path to the created sorting file for sorting ffindex databases.
    :rtype: pathlib.Path
    """
    if sorting_file is None:
        sorting_file = db_dir.joinpath("sorting.dat")

    # The parameter set is supplied by the HHsuite3 documentation
    sort_command = f"sort -k3 -n -r '{cs219_ffindex}' "
    cut_command = "cut -f1"

    split_sort_command = shlex.split(sort_command)
    split_cut_command = shlex.split(cut_command)

    # Prepares to pipe the output of the sort command into cut
    # Cut grabs the first field and writes it out
    sort_process = Popen(args=split_sort_command, stdout=PIPE)
    cut_process = Popen(args=split_cut_command, stdout=PIPE,
                        stdin=sort_process.stdout)
    sort_process.stdout.close()

    # Writes the output of the cut command into file, after decoding
    with sorting_file.open(mode="w") as filehandle:
        out, err = cut_process.communicate()
        filehandle.write((out).decode("utf-8"))

    return sorting_file


def sort_ffindex(db_dir, sorting_file, ffindex_tuple, stdout=None):
    """Function to sort a ffindex database with a sorting key.

    :param db_dir: Path to a directory where the database will be created.
    :type db_dir: pathlib.Path
    :param sorting_file: Path to a sorting file key for ffindex databases.
    :type sorting_file: pathlib.Path
    :param ffindex_tuple: A tuple containing path objects to a ffindex db.
    :type ffindex_tuple: tuple(pathlib.Path, pathlib.Path)
    """
    if stdout is None:
        stdout = DEVNULL

    # Defines the path for the ordered ffindex database
    ord_ffdata = db_dir.joinpath("".join([ffindex_tuple[0].stem,
                                          "ordered.ffdata"]))
    ord_ffindex = db_dir.joinpath("".join([ffindex_tuple[1].stem,
                                           "ordered.ffindex"]))

    command = (f"ffindex_order '{sorting_file}' "
               f"'{ffindex_tuple[0]}' '{ffindex_tuple[1]}' "
               f"'{ord_ffdata}' '{ord_ffindex}'")

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    # Overwrites the unsorted ffindex database files with sorted files
    ord_ffdata.replace(ffindex_tuple[0])
    ord_ffindex.replace(ffindex_tuple[1])


# Placeholder function
def verify_hhsuite_database(db_dir, db_name, cores=1, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    return True


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def create_version_file(filepath, versions_dict):
    with filepath.open(mode="w") as filehandle:
        json.dump(versions_dict, filehandle)
