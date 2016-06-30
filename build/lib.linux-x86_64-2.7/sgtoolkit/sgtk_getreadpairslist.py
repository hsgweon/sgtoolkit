#!/usr/bin/env python

""" 
Produces a list of filenames for each read pairs and their library names from demultiplxed illumina fastq.gz files
"""

import sys, os, subprocess, ConfigParser, shutil
from time import strftime

__author__ = "Hyun Soon Gweon"
__copyright__ = "Copyright 2015, The PIPITS Project"
__credits__ = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__ = "GPL"
__maintainer__ = "Hyun Soon Gweon"
__email__ = "hyugwe@ceh.ac.uk"

HEADER = "\033[95m"
BLUE = "\033[94m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
CYAN = "\033[96m"
RED = "\033[91m"
ENDC = "\033[0m"

def make_read_pairs_list(options):

    print(GREEN + "Generating a read-pair list file from the input directory..." + ENDC)

    samples = []
    fastqs = []
    for file in os.listdir(options.inputdir):
        if file.endswith(".fastq.gz") or file.endswith(".fastq.bz2") or file.endswith(".fastq"):
            samples.append(file.split("_")[0])
            samples = list(set(samples))
            fastqs.append(file)
    
    # Check files
    offendingSampleIDs = []
    for sample in samples:
        if len([s for s in fastqs if sample in s]) != 2:
            offendingSampleIDs.append(sample)

    if len(offendingSampleIDs) != 0:
        print(RED + "N.B. There are duplicate/missing pair(s) in the Illumina sequences. Correct your filenames before continuing." + ENDC)
        print(RED + "N.B. Offending sample ID(s): " + ", ".join(offendingSampleIDs) + ENDC)

    
    # All passed
    if len(fastqs) % 2 != 0:
        print(RED + "There are missing pair(s) in the Illumina sequences. Check your files and labelling" + ENDC)
        exit(1)

    fastqs_l = []
    fastqs_f = []
    fastqs_r = []

    coin = True
    for fastq in sorted(fastqs):
        if coin == True:
            fastqs_f.append(fastq)
        else:
            fastqs_r.append(fastq)
        coin = not coin

    for i in range(len(fastqs_f)):
        if fastqs_f[i].split("_")[0] != fastqs_r[i].split("_")[0]:
            logger.error("Problem with labelling the files.")
            exit(1)
        fastqs_l.append(fastqs_f[i].split("_")[0])

    if not options.output:
        outfile_name = "readpairslist.txt"
    else:
        outfile_name = options.output

    outfile_fastqslist = open(outfile_name, "w")
    outfile_fastqslist.write("# Lines beginning with \"#\" is ignored. \n")
    outfile_fastqslist.write("# SampleID\tFilename for forward reads\tFilename for reverse reads\n")
    count = 1
    for i in range(len(fastqs_f)):
        label = ""
        if options.label_add_c:
            label = fastqs_l[i] + options.label_add_c
        elif options.label_add_reindex_c:
            label = options.label_add_reindex_c + str(count).zfill(3)
        else:
            label = fastqs_l[i]
        count += 1
        outfile_fastqslist.write(label + "\t" + fastqs_f[i] + "\t" + fastqs_r[i] + "\n")
    outfile_fastqslist.close()

    print("Done. \"" + outfile_name + "\" created.")


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("makes a read_pairs_list.")
    parser.add_argument(
        "-i",
        action = "store",
        dest = "inputdir",
        metavar = "<DIR>",
        help = "[REQUIRED] Directory with your raw sequences in gzipped FASTQ",
        required = True)
    parser.add_argument(
        "-o",
        action = "store",
        dest = "output",
        metavar = "<FILE>",
        help = "Name of output list file.",
        required = False)
    parser.add_argument(
        "--label-add-c",
        action = "store",
        dest = "label_add_c",
        metavar = "<TXT>",
        help = "Add a label to each sample ids in the output file. N.B. \"_\" is not allowed",
        required = False)
    parser.add_argument(
        "--label-reindex-c",
        action = "store",
        dest = "label_add_reindex_c",
        metavar = "<TXT>",
        help = "Rename samples with the given label. It will automatically add 001, 002 etc. at the end of each name. N.B. \"_\" is not allowed",
        required = False)
    options = parser.parse_args()

    if options.label_add_c:
        if options.label_add_c.find("_") != -1:
            print("Error: \"_\" is not allowed in the sample id")
            exit(1)
    
    make_read_pairs_list(options)

