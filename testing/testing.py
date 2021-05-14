"""Testing module for Prophicient with metrics and statistics."""

from Bio import SeqIO
from pathlib import Path

import csv
import argparse
import math


# Global variables
CSV_HEADERS = ["strain", "software", "prophage", "ends"]    # headers for csv
GENOME_LENGTHS = {}     # stores lengths of the genomes


# Will convert the file to a biopython object
def open_file(filename, store_length=False):
    """Convert the FASTA formatted file into a biopython object.

    :param filename: name of the file
    :type filename: Path
    :param store_length: True if the lengths need to be stored
    :type store_length: bool
    :return: SeqRecord object generated from the FASTA file
    :rtype: Bio.SeqRecord.SeqRecord
    """
    # iterate through all the contigs
    # this makes sure even extrachromosomal prophages are seen
    contigs = [record for record in SeqIO.parse(filename, "fasta")]

    length = 0

    # length of the genome = sum of lengths of every record
    for record in contigs:
        length += len(record)

    # store lengths of the genomes in the dictionary
    if store_length:
        GENOME_LENGTHS[filename.stem] = length

    return contigs


def find_sub(parent, child_contigs):
    """Locate the prophage region as a subsequence of the bacterial genome.

    :param parent: Parent sequence
    :type parent: Bio.SeqRecord.SeqRecord
    :param other: Subsequence being located in the parent
    :type other: Bio.SeqRecord.SeqRecord
    :return: left end of the subsequence
    :rtype: int
    """
    # keep track of which contig it is in

    for record in child_contigs:
        # Find the beginning of the subsequence
        left_coor = parent.seq.find(record.seq)
        # Checks reverse complement
        if left_coor == -1:
            reverse = record.reverse_complement()
            left_coor = parent.seq.find(reverse.seq)
            # if nothing is returned an error exists
            if left_coor == -1:
                return None
    return left_coor


def stats(manual_data, testing_data):
    """Return True Positive, False Positive, True Negative, False Negative.

    :param manual_coordinates: Manual prophage names w/ Coordinates
    :type manual_coordinates: dict
    :param software_coordinates: Software prophages names w/ Coordinates
    :type software_coordinates: dict
    :return: Dictionary of TP, TF, FP and FN values
    :rtype: dict
    """
    statistics = {}
    statistics["TRUE_POSITIVE"] = 0
    statistics["FALSE_POSITIVE"] = 0
    statistics["FALSE_NEGATIVE"] = 0
    manual_set = set()
    soft_set = set()

    for ends in testing_data.values():  # iterates over the ends
        ends = eval(ends)
        for coor in range(int(ends[0]), int(ends[1])):
            soft_set.add(coor)
    for ends in manual_data.values():  # iterates over the ends for manual
        ends = eval(ends)
        for coor in range(int(ends[0]), int(ends[1])):
            manual_set.add(coor)

    statistics["TRUE_POSITIVE"] += len(soft_set.intersection(manual_set))
    statistics["FALSE_POSITIVE"] += len(soft_set.difference(manual_set))
    statistics["FALSE_NEGATIVE"] += len(manual_set.difference(soft_set))
    return statistics


def collect_stats(filepath):
    """Compile raw bp values of TP,FP,TN,FN.

    :param filepath: Filepath of the csv file with data entries
    :type filepath: Path
    :return: Dictionary containing TP, TN, FP, FN
    :rtype: dict
    """
    true_positive = 0  # the bp of true positives
    true_negative = 0  # the bp of true negatives
    false_positive = 0  # the bp of false positives
    false_negative = 0  # the bp of
    length_manual = 0
    length_test = 0
    manual_ends = {}  # {name: (start, stop)}
    testing_ends = {}
    statistics = {}

    statistics["TRUE_POSITIVE"] = 0
    statistics["FALSE_POSITIVE"] = 0
    statistics["FALSE_NEGATIVE"] = 0
    statistics["TRUE_NEGATIVE"] = 0

    with filepath.open(mode="r") as filehandle:

        csv_reader = csv.reader(filehandle, delimiter=",", quotechar='"')

        for strain in GENOME_LENGTHS:
            for row in csv_reader:
                # for manual data
                if row[1] == "manual" and row[0] == strain:
                    manual_ends[row[2]] = row[3]
                    # convert to tuple
                    ends_tuple = eval(manual_ends.get(row[2]))
                    # check forward or reverse
                    if ends_tuple[1] - ends_tuple[0] < 0:
                        length_manual += ends_tuple[0] - ends_tuple[1]
                    else:
                        length_manual += ends_tuple[1] - ends_tuple[0]
                # for software data - same logic as above
                if row[1] == filepath.stem.split(
                        "_")[0] and row[0] == strain:
                    testing_ends[row[2]] = row[3]
                    # convert to tuple
                    test_ends = eval(testing_ends.get(row[2]))
                    if test_ends[1] - test_ends[0] < 0:
                        length_test += test_ends[0] - test_ends[1]
                    else:
                        length_test += test_ends[1] - test_ends[0]

            per_phage = stats(manual_ends, testing_ends)
            true_positive += per_phage.get("TRUE_POSITIVE")
            false_positive += per_phage.get("FALSE_POSITIVE")
            false_negative += per_phage.get("FALSE_NEGATIVE")

            true_negative += GENOME_LENGTHS[strain] - \
                length_manual - false_positive
    statistics["TRUE_POSITIVE"] += true_positive
    statistics["FALSE_POSITIVE"] += false_positive
    statistics["FALSE_NEGATIVE"] += false_negative
    statistics["TRUE_NEGATIVE"] += true_negative
    return statistics


def metrics(stats):
    """Calculate the metric data.

    :param stats: statistics data with tp, tn, fn, fp values
    :type stats: dict
    :return: Dictionary containing sensitivty, positive predictive value,
             accuracy and Matthews Correlation Coefficient (MCC)
    :rtype: dict
    """
    # collect data
    tp = stats["TRUE_POSITIVE"]
    tn = stats["TRUE_NEGATIVE"]
    fp = stats["FALSE_POSITIVE"]
    fn = stats["FALSE_NEGATIVE"]

    # calculate metric values
    sn = (tp/(tp+fn))*100
    ppv = (tp/(tp+fp))*100
    acc = ((tp+tn)/(tp+tn+fp+fn))*100
    mcc = ((tn*tp) - (fn*fp))/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

    # return metrics as dictionary
    return {"sensitivity": round(sn, 3),
            "ppv": round(ppv, 3),
            "accuracy": round(acc, 3),
            "mcc": round(mcc, 3)}


def get_child_ends(prophage_filepath, parent_contigs):
    """Return the ends of the prophage region as a tuple.

    :param prophage_filepath: Filepath to the hypothetical prophage FASTA
    :type prophage_filepath: Path
    :param parent_seq: Sequence of the parent genome
    :return: Tuple of prophage ends
    :rtype: tuple
    """
    # read all the prophages in that folder

    child_contigs = open_file(prophage_filepath)
    child_tuples = []

    for child_seq in child_contigs:
        for parent_seq in parent_contigs:
            # find the start of the prophage regions - contigs
            start = find_sub(parent_seq, child_contigs)

            # not found
            if start is None:
                continue

            # tuple with ends (stop = start+length)
            child_tuples.append((start, start+len(child_seq)))

    return child_tuples


def comparison(dir, manual_name, software_name):
    """Compare the data.

    :param dir: Path to working directory
    :type dir: Path
    :param manual_name: Name of directory with manually annotated genomes
    :type manual_name: str
    :param software_name: Name of software being tested
    :type software_name: str
    """
    # user enters the name of the folder with manually annotated
    # and the software name
    path_to_manual = dir/"prophages"/manual_name
    path_to_software = dir/"prophages"/software_name

    genome_path = dir/"genomes"

    csv_filename = software_name + "_testing_results.csv"

    csv_path = dir/csv_filename

    with open(csv_path, mode="w") as csv_file:

        writer = csv.DictWriter(csv_file, fieldnames=CSV_HEADERS)
        writer.writeheader()

        for genome in sorted(genome_path.iterdir()):
            strain_name = genome.stem
            ref_seq = open_file(genome, True)     # reference genome sequence
            manual_strain_path = path_to_manual/strain_name
            software_strain_path = path_to_software/strain_name

            # data to csv is written as dictionaries with headers as keys
            for prophage in sorted(manual_strain_path.iterdir()):
                manual_ends = get_child_ends(prophage, ref_seq)
                if manual_ends is None:
                    # this should not happen
                    print("ERROR: Prophage not found: ", prophage)
                    dict = {"strain": strain_name,
                            "software": "manual",
                            "prophage": prophage.stem,
                            "ends": "NONE"}
                    writer.writerow(dict)
                else:
                    # data written for each manually annotated genome
                    for ends_tuple in manual_ends:
                        dict = {"strain": strain_name,
                                "software": "manual",
                                "prophage": prophage.stem,
                                "ends": ends_tuple}
                        writer.writerow(dict)

            for prophage in sorted(software_strain_path.iterdir()):
                # check that no hidden files are being read in
                # only read in fasta files
                if prophage.suffix == ".fasta":
                    software_ends = get_child_ends(prophage, ref_seq)
                else:
                    continue

                # same logic as above for software data
                if software_ends is None:
                    # if the prophage found by the software not found in host
                    # this should not happen
                    dict = {"strain": strain_name,
                            "software": software_name,
                            "prophage": prophage.stem,
                            "ends": "NONE"}
                    writer.writerow(dict)
                else:
                    # write software data to the csv
                    for ends_tuple in software_ends:
                        dict = {"strain": strain_name,
                                "software": software_name,
                                "prophage": prophage.stem,
                                "ends": ends_tuple}
                        writer.writerow(dict)
        csv_file.close()    # close csv file

        stats = collect_stats(csv_path)     # collect the statistics
        metric_data = metrics(stats)        # metrics
        print_data(stats, metric_data)      # print all the relevant data


def print_data(stats, metrics):
    """Print the stats and metric data.

    :param stats: Statistics for the software
    :type stats: dict
    :param metrics: Metrics for the software
    :type metrics: dict
    """
    print("--------------------------------------------------------------")
    print("TRUE POSITIVE\tFALSE POSITIVE\tTRUE NEGATIVE\tFALSE NEGATIVE")
    print("--------------------------------------------------------------")
    print(f"{stats['TRUE_POSITIVE']}\t\t{stats['FALSE_POSITIVE']}\t",
          f"\t{stats['TRUE_NEGATIVE']}\t{stats['FALSE_NEGATIVE']}")

    print("\n\n-----------------------------------------------------------")
    print("Sensitivity\tPPV\t\tAccuracy\tMCC")
    print("-----------------------------------------------------------")
    print(f"{metrics['sensitivity']}\t",
          f"\t{metrics['ppv']}\t\t{metrics['accuracy']}\t\t{metrics['mcc']}")


def get_args():
    """Parse the argument list.

    :return: Working directory where the sequences are located
    :rtype: Path
    """
    parser = argparse.ArgumentParser(description=__doc__)
    dir_help = "Path to working directory"
    manual_help = "Name of directory with manually annotated prophages"
    software_help = "Name of the software being tested"

    parser.add_argument("dir", type=Path, help=dir_help)
    parser.add_argument("manual", type=str, help=manual_help)
    parser.add_argument("software", type=str, help=software_help)
    args = parser.parse_args()

    return args


def main():
    """Run the testing module."""
    # 3 args: reference path, manual path, software path
    args = get_args()

    print(f"\n{args.software}")     # print software name

    comparison(args.dir, args.manual, args.software)    # run comparison


if __name__ == '__main__':
    main()
