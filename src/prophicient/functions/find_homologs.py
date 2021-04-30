from prophicient.classes import hhresult_parsing
from prophicient.functions import multiprocess, wrapper_basic


DEFAULTS = {"expect_cutoff": 0.0001, "probability": 90}


def find_homologs(trans_dir, output_dir, database_path, cores=1,
                  expect_cutoff=DEFAULTS["expect_cutoff"],
                  probability=DEFAULTS["probability"],
                  sort_attr="score", verbose=False):
    """Searches an HHsuite database for every file in the input gene
    translations directory and returns the translation ID and
    the ID of the HMM profile of it's best alignment.

    :param trans_dir: Directory of files to perform an HHsearch on.
    :type trans_dir: pathlib.Path
    :param output_dir: Directorry to place HHsearch alignment result files.
    :type output_dir: pathlib.Path
    :param database: Path to the HHsuite database.
    :type database: pathlib.Path
    :param cores: Number of HHsearch processes to create.
    :type cores: int
    :param expect_cutoff: E-value cutoff to limit HHsearch alignment matches.
    :type expect_cutoff: float
    :param probability: Probability cutoff to limit HHsearch alignment matches.
    :type probability: float
    :param sort_attr: Sort key of HHsearch alignment matches.
    :type sort_attr: str
    :param verbose: Toggle verbose print statements of the function.
    :type verbose: bool
    """
    # Create a list of arguments to pass to the HHsearch wrapper function
    work_items = create_job_queue(trans_dir, output_dir, database_path,
                                  expect_cutoff)

    # Parallelize HHsearch calls
    multiprocess.parallelize(work_items, cores, wrapper_basic.hhsearch,
                             verbose=verbose)

    hhresult_objects = []
    for hhresult_file in output_dir.iterdir():
        # Apply probability cutoff to HHsearch results
        # If the result passes the given threshold, parse it into an object
        hhresult = test_homology(hhresult_file, probability)
        if hhresult is not None:
            hhresult_objects.append(hhresult)

    for hhresult_object in hhresult_objects:
        # Sort the HHresult matches with the given sort key
        hhresult_object.matches.sort(key=lambda x: int(getattr(x, sort_attr)),
                                     reverse=True)

    # Take the first HHresult match from the sorted matches
    homologs = [(hhresult.query_id, hhresult.matches[0].target_id)
                for hhresult in hhresult_objects]

    return homologs


def create_job_queue(input_dir, output_dir, database, cutoff):
    """Creates a job queue for hhsearch items from an input and output
    directory.

    :param input_dir: Directory of files to perform an HHsearch on.
    :type input_dir: pathlib.Path
    :param output_dir: Directory to place HHsearch alignment result files.
    :type output_dir: pathlib.Path
    :param database: Path to the HHsuite database.
    :type database: pathlib.Path
    :param cutoff: E-value cutoff to limit HHsearch alignment matches.
    :type cutoff: float
    :return: Work items to parallelize HHsearch calls.
    :rtype: list[tuple]
    """
    jobs = []

    for filepath in input_dir.iterdir():
        if filepath.suffix == ".fasta":
            jobs.append((filepath, output_dir, database, cutoff))

    return jobs


def test_homology(hhresult_file, probability):
    """Tests whether the most probable hit of an HMM profile alignment meets
    the given probability alignment.

    :param hhresult_file: Filepath to a hhsearch alignment result.
    :type hhresult_file: pathlib.Path
    :param probability: Probability threshold.
    :type probability: pathlib.Path
    :return: A valid HHResult object if the given alignment meets the threshold
    :rtype: prophicient.classes.hhresult_parsing.HHResult
    """
    # Parse the HHSearch alignment result file into an object
    hhresult = hhresult_parsing.HHResult(hhresult_file)
    hhresult.parse_result()

    if not hhresult.matches:
        return None

    # Sort the hhresult matches by probability, and compare the highest
    # probability alignment against the threshold
    hhresult.matches.sort(key=lambda x: x.probability, reverse=True)
    if float(hhresult.matches[0].probability) < probability:
        return None

    return hhresult
