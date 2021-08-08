"""
Functions to parallelize of processing of a list of inputs.
Adapted from https://docs.python.org/3/library/multiprocessing.html
"""
import multiprocessing as mp

from depht.classes.progress import Progress

# Make sure new processes are forked, not spawned
mp.set_start_method("fork")

# Assume everyone has hyper-threading; only use physical cores
LOGICAL_CORES = mp.cpu_count()
PHYSICAL_CORES = LOGICAL_CORES // 2


def show_progress(current, end, width=50):
    """Updating progress bar that prints to the command line in-line.

    :param current: current value (0-n)
    :param end: end value (n)
    :param width:
        character width for the progress bar (default 50);
        ideally 100 divided by width should naturally be an integer
    :return:
        progress - the integer percent completion calculated
        as the ratio of current to end multiplied by 100
    """
    progress = Progress(current, end, width)
    print(f"\r{progress}", end="")
    return progress


def parallelize(inputs, num_processors, task, verbose=False):
    """
    Parallelize some task on an input list across the specified number
    of processors
    :param inputs: list of inputs
    :param num_processors: number of processor cores to use
    :param task: name of the function to run
    :param verbose: updating progress bar output?
    :return: results
    """
    results = []

    # Don't do any work if there are no inputs
    if len(inputs) == 0:
        return results

    num_processors = count_processors(inputs, num_processors)

    tasks = []
    for item in inputs:
        if not isinstance(item, tuple):
            item = (item,)
        tasks.append((task, item))

    # Start working on the jobs
    results = start_processes(tasks, num_processors, verbose)

    return results


def count_processors(inputs, num_processors):
    """
    Programmatically determines whether the specified num_processors is
    appropriate. There's no need to use more processors than there are
    inputs, and it's impossible to use fewer than 1 processor or more
    than exist on the machine running the code.
    :param inputs: list of inputs
    :param num_processors: specified number of processors
    :return: num_processors (optimized)
    """
    if num_processors < 1 or num_processors > LOGICAL_CORES:
        print(f"Invalid number of CPUs specified ({num_processors})")
        num_processors = min([mp.cpu_count(), len(inputs)])
        print(f"Using {num_processors} CPUs...")

    return num_processors


def worker(input_queue, output_queue):
    for func, args in iter(input_queue.get, 'STOP'):
        result = func(*args)
        output_queue.put(result)
    return


def start_processes(inputs, num_processors, verbose):
    """
    Creates input and output queues, and runs the jobs
    :param inputs: jobs to run
    :param num_processors: optimized number of processors
    :param verbose: updating progress bar output?
    :return: results
    """
    job_queue = mp.Queue()
    done_queue = mp.Queue()

    # Counter so we know how many tasks we have in all (input + show_progress tasks)
    tasks = 0

    if verbose is True:
        interval = max([1, len(inputs)//100])

        # Put inputs into job queue
        for i in range(len(inputs)):
            tasks += 1
            if i % interval == 0:
                job_queue.put((show_progress, (i, len(inputs))))
                tasks += 1
            job_queue.put(inputs[i])
        tasks += 1
        job_queue.put((show_progress, (len(inputs), len(inputs))))
    else:
        for i in range(len(inputs)):
            tasks += 1
            job_queue.put(inputs[i])

    # Put a bunch of 'STOP' signals at the end of the queue
    for i in range(num_processors):
        job_queue.put('STOP')

    # Start up workers
    worker_pool = []
    for i in range(num_processors):
        worker_n = mp.Process(target=worker, args=(job_queue, done_queue))
        worker_n.start()
        worker_pool.append(worker_n)

    # Grab results from done queue
    results = []

    # Remove non-Progress results
    for i in range(tasks):
        result = done_queue.get()
        if not isinstance(result, Progress):
            results.append(result)

    [worker_n.join() for worker_n in worker_pool]

    # Leave the progress bar line
    if verbose is True:
        print("\n")

    return results
