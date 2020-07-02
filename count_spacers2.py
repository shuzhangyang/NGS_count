'''
Count sgRNAs from NGS fastq data 
Developed from Joung 2017 Feng Zhang etal. Supplementary Data 3: count_spacers.py
'''

from Bio import SeqIO
import csv
from collections import OrderedDict
import numpy as np
import argparse
import os, glob
from timeit import default_timer as timer

KEY_REGION_START = 30 #start index of key region
KEY_REGION_END = 55 #end index of key region
KEY0 = "CGAAACACC" #identifies sequence before guide to determine guide position

def count_spacers(lib_file, fastq_file, output_file, guide_g): 
    """
    creates a dictionary with guide counts from fastq_file, writes to output_file
    fastq_file: forward read fastq file
    output_file: csv file to write guide dictionary to
    dictionary: guide sequence as key, guide count as entry
    """

    num_reads = 0 #total number of reads processed
    perfect_matches = 0 # guides with perfect match to library
    non_perfect_matches = 0 #number of guides without a perfect match to the library
    key_not_found = 0 #count of reads where key was not found

    # add 'G' to key sequence if included in library
    if guide_g:
#        global KEY
        KEY = KEY0+"G"

    # open library sequences and initiate dictionary of read counts for each guide
    try:
        with open(lib_file) as infile: 
            reader = csv.reader(infile)
            dictionary = {rows[0]:0 for rows in reader}
    except:
        print('could not open', lib_file)
      
    # open fastq file
    try:
        handle = open(fastq_file)
    except:
        print("could not find fastq file")
        return

    # process reads in fastq file
    readiter = SeqIO.parse(handle, "fastq")
    for record in readiter: #contains the seq and Qscore etc.
        num_reads += 1
        read_sequence = str.upper(str(record.seq))
        key_region = read_sequence[KEY_REGION_START:KEY_REGION_END]
        key_index = key_region.find(KEY)
        if key_index >= 0:
            start_index = key_index + KEY_REGION_START + len(KEY)
            guide = read_sequence[start_index:(start_index + 20)]
            if guide in dictionary:
                dictionary[guide] += 1
                perfect_matches += 1
            else:
                non_perfect_matches += 1
        else:
            key_not_found += 1

    # create ordered dictionary with guides and respective counts and output as a csv file                      
    dict_sorted = OrderedDict(sorted(dictionary.items(), key=lambda t: t[0]))
    with open(output_file, 'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',')
        for guide in dict_sorted:
            count = dict_sorted[guide]
            mywriter.writerow([guide,count])

    # percentage of guides that matched perfectly
    percent_matched = round(perfect_matches/float(perfect_matches + non_perfect_matches) * 100, 1)
    # percentage of undetected guides with no read counts
    guides_with_reads = np.count_nonzero(list(dictionary.values()))
    guides_no_reads = len(list(dictionary.values())) - guides_with_reads
    percent_no_reads = round(guides_no_reads/float(len(list(dictionary.values()))) * 100, 1)
    # skew ratio of top 10% to bottom 10% of guide counts
    top_10 = np.percentile(list(dictionary.values()), 90)
    bottom_10 = np.percentile(list(dictionary.values()), 10)
    if top_10 != 0 and bottom_10 != 0:
        skew_ratio = top_10/bottom_10
    else:
        skew_ratio = 'Not enough perfect matches to determine skew ratio'

    # write analysis statistics to output_file_statistics.csv
    with open(os.path.splitext(fastq_file)[0]+'_stat.csv', 'w') as stat:

                statwriter = csv.writer(stat, delimiter=',')
                statwriter.writerow(['perfect_matches', perfect_matches])
                statwriter.writerow(['non_perfect_matches', non_perfect_matches])
                statwriter.writerow(['key_not_found', key_not_found])
                statwriter.writerow(['num_reads', num_reads])
                statwriter.writerow(['percent_matched', percent_matched])
                statwriter.writerow(['percent_no_reads', percent_no_reads])
                statwriter.writerow(['skew_ratio', skew_ratio])

    handle.close()           
    return 
    

if __name__ == '__main__':
    start = timer()
    
    parser = argparse.ArgumentParser(description='Analyze sequencing data for sgRNA library distribution')
    parser.add_argument('-f', '--fastq', type=str, dest='fastq_file',
                        help='fastq file name', default='NGS.fastq')
    parser.add_argument('-o', '--output', type=str, dest='output_file',
                        help='output file name', default='library_count.csv')
    parser.add_argument('-i', '--input', type=str, dest='lib_file',
                        help='input file name', default='library_sequences.csv')
    parser.add_argument('-no-g', dest='guide_g', help='presence of guanine before spacer', action='store_false')
    parser.set_defaults(guide_g=True)
    args = parser.parse_args()
    
    #loop through fastq files
    files = glob.glob("*.fastq")
    for file in files:
        count_spacers('lib_brie_sgRNA.csv', file, os.path.splitext(file)[0]+'_Count.csv', True)

    print('Total time = '+str(timer()-start))
    