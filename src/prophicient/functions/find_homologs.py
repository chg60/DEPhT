from prophicient.classes import hhresult_parsing
from prophicient.functions import multiprocess, wrapper_basic


DEFAULTS = {"expect_cutoff": 0.0001, "probability": 90}


def find_homologs(trans_dir, output_dir, database_path, cores=1,
                  expect_cutoff=DEFAULTS["expect_cutoff"],
                  probability=DEFAULTS["probability"],
                  sort_attr="score", verbose=False):
    work_items = create_job_queue(trans_dir, output_dir, database_path,
                                  expect_cutoff)

    multiprocess.parallelize(work_items, cores, wrapper_basic.hhsearch,
                             verbose=verbose)

    hhresult_objects = []
    for hhresult_file in output_dir.iterdir():
        hhresult = test_homology(hhresult_file, probability)
        if hhresult is not None:
            hhresult_objects.append(hhresult)

    hhresult_objects.sort(key=lambda x: getattr(x.matches[0], sort_attr),
                          reverse=True)

    homologs = [(hhresult.query_id, hhresult.matches[0].target_id)
                for hhresult in hhresult_objects]

    return homologs


def create_job_queue(input_dir, output_dir, database, cutoff):
    jobs = []

    for filepath in input_dir.iterdir():
        if filepath.suffix == ".fasta":
            jobs.append((filepath, output_dir, database, cutoff))

    return jobs


def test_homology(hhresult_file, probability):
    hhresult = hhresult_parsing.HHResult(hhresult_file)
    hhresult.parse_result()

    if not hhresult.matches:
        return None

    hhresult.matches.sort(key=lambda x: x.probability, reverse=True)
    if float(hhresult.matches[0].probability) < probability:
        return None

    return hhresult
