from Bio.Data.CodonTable import TranslationError
from Bio import SeqIO
from Bio.Seq import Seq
import sys

GENE_CUTOFF = 0.8

def read_file(filename):
        seq_record=SeqIO.read(filename, "gb")
        return seq_record


def find_stop(seq_obj, prod_obj):
    stop_dictionary={}
    start_and_stop_dictionary={}
    start_dictionary={}
    total_count=0
    true_count=0
    prod_count=0

    for feature in seq_obj.features:

        if feature.type  not in ("CDS","tRNA", "tmRNA"):
            continue
        total_count+=1
        start=feature.location.nofuzzy_start
        stop=feature.location.nofuzzy_end

        # manual file
        stop_dictionary[stop]=False
        start_and_stop_dictionary[stop]=start;
        start_dictionary[start]=0;

    for feature in prod_obj.features:
        if feature.type  not in ("CDS","tRNA", "tmRNA"):
            continue
        next_start=feature.location.nofuzzy_start
        next_stop=feature.location.nofuzzy_end

        #if feature.strand==-1:
            #temp=next_start
            #next_start=next_stop
            #next_stop=temp

        if stop_dictionary.get(next_stop)!=None:
            true_count+=1 # true positive
            new_key= start_and_stop_dictionary.get(next_stop)
            start_dictionary[new_key]=next_start
        else:
            prod_count+=1

    print(true_count) # true_count = number of matching stops
    print(total_count) #total_count = all the genes

    test_output = []

    # comparison of number of genes called in both files
    if (true_count/total_count)>GENE_CUTOFF: # rough estimate - 80%
        test_output.append(True) # if prodigal returns 80% of genes
    else:
        test_output.append(False) # otherwise False

    true_positive = true_count
    false_positive = prod_count # in prodigal but not in manual
    false_negative = total_count - true_count # in manual file but not in prodigal output


    test_output.append(true_positive)
    test_output.append(false_positive)
    test_output.append(false_negative)

    return test_output

def sensitivity(true_positive, false_negative):
    sn = (true_positive)/(true_positive+false_negative)
    sn *= 100
    return sn

def positive_predictive_value(true_positive, false_positive):
    ppv = (true_positive)/(true_positive+false_positive)
    ppv *= 100
    return ppv

    # feed in tp, fp, fn
    # calculate sn, ppv, accuracy

def print_output(data):

    print("> ", GENE_CUTOFF, " genes detected = ", data[0])
    print("sesitivity = ", sensitivity(data[1], data[3]))
    print("positive predictive value = ", positive_predictive_value(data[1], data[2]))



def main():
    prod_file=sys.argv[1] # prodigal
    og_file=sys.argv[2] # manual file
    prod_obj=read_file(prod_file)
    og_obj=read_file(og_file)
    data=find_stop(og_obj, prod_obj)

    print_output(data)

    # print(results[0])



if __name__ == '__main__':
    main()
