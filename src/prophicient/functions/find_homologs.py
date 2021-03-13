from prophicient.classes import hhresult_parsing
from prophicient.functions import multiprocesses, wrapper_basic


DEFAULTS = {"expect_cutoff": 0.0001, "probability": 0.9}


def find_homologs(trans_dir, output_dir, database_path, cores=1,
                  expect_cutoff=DEFAULTS["expect_cutoff"],
                  probability=DEFAULTS["probability"],
                  verbose=False):
    work_items = []
    for input_file in trans_dir.iter_dir():
        work_items.append((input_file,))

    multiprocesses.parallelize(work_items, cores, wrapper_basic.hhsearch,
                               verbose=verbose)

    hhresult_objects = []
    for hhresult_file in output_dir.iterdir():
        hhresult = test_homology(hhresult_file, probability)
        if hhresult is not None:
            hhresult_objects.append(hhresult)

    return hhresult_objects


def test_homology(hhresult_file, probability):
    hhresult = hhresult_parsing.HHResult(hhresult_file)
    hhresult.parse_result()

    if not hhresult.matches:
        return None

    hhresult.matches.sort(key=lambda x: x.probability, reverse=True)
    if hhresult.matches[0].probability < probability:
        return None

    return hhresult
