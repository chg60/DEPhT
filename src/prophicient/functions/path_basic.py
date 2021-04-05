import sys
from pathlib import Path


def make_new_dir(output_dir, new_dir, attempt=1, mkdir=True):
    """Make a new directory.

    Checks to verify the new directory name is valid and does not
    already exist. If it already exists, it attempts to extend
    the name with an integer suffix.

    :param output_dir:
        Full path to the directory where the new directory will be created.
    :type output_dir: Path
    :param new_dir: Name of the new directory to be created.
    :type new_dir: Path
    :param attempt: Number of attempts to create the directory.
    :type attempt: int
    :returns:
        If successful, the full path of the created directory.
        If unsuccessful, None.
    :rtype: Path, None
    """
    valid = False
    count = 0
    while (not valid and count < attempt):
        if count > 0:
            new_dir_mod = new_dir.stem + "_" + str(count)
            new_dir_mod = Path(new_dir_mod)
        else:
            new_dir_mod = new_dir
        new_path = Path(output_dir, new_dir_mod)
        if not new_path.is_dir():
            valid = True
            if mkdir:
                new_path.mkdir()
        count += 1
    if not valid:
        return None
    else:
        return new_path


def create_working_dir(working_path, dump=False, force=False):
    if not dump:
        if not force:
            if working_path.is_dir():
                print("COWARDLY ABORTING PIPELINE: "
                      f"Directory '{working_path}' already exists.")
                sys.exit(1)

    working_path.mkdir(parents=True, exist_ok=True)


def create_working_path(folder_path, folder_name, dump=False, force=False,
                        attempt=50):
    if folder_path is None:
        working_path = create_default_path(folder_name, force=force,
                                           attempt=attempt)
    else:
        working_path = folder_path.joinpath(folder_name)

    if dump:
        working_path = working_path.parent

    return working_path


def create_default_path(name, force=False, attempt=50):
    default_path = Path.cwd().joinpath(name)
    if not force:
        default_path = make_new_dir(Path.cwd(), default_path,
                                    attempt=attempt, mkdir=False)

    return default_path


def convert_file_path(path_string):
    """Function to convert argparse input to a working file path.

    :param path: A string to be converted into a Path object.
    :type path: str
    :returns: A Path object converted from the inputted string.
    :rtype: pathlib.Path
    """
    path = Path(path_string)
    path = path.expanduser()
    path = path.resolve()

    if not path.is_file():
        print(f"The file {path} does not exist.")
        sys.exit(1)

    return path
