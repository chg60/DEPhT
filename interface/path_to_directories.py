"""Create a list of files for the GUI."""

from pathlib import Path


def extract_files(linear_dir, circular_dir):
    """Extract all the paths from the given directories.

    :param linear_dir: name of folder w/ annotated gb files
    :type linear_dir: str
    :param circular_dir: name of folder w/ prophage features
    :type linear_dir: str
    :return: list of all the paths needed for display
    :rtype: list
    """
    dir_name = Path.cwd()
    lin_path_name = dir_name / linear_dir
    cir_path_name = dir_name / circular_dir

    paths = []
    path_circular = []
    path_linear = []

    for file in sorted(lin_path_name.iterdir()):
        # NAME = str(file.name)
        if file.suffix == ".gb":
            stem = str(file)
            path_linear.append(stem)
    paths.append(path_linear)
    for file in cir_path_name.iterdir():
        if file.suffix == ".gb":
            stem = str(file)
            path_circular.append(stem)
    paths.append(path_circular)
    # print(paths)

    return paths


def main():
    """Run the function."""
    extract_files("prophages", "output")


main()
