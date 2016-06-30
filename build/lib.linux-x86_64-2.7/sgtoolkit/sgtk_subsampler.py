#!/usr/bin/env python

####################
# Argument Options #
####################

import argparse
parser = argparse.ArgumentParser(description = "Randomly subsamples sequences")
parser.add_argument("-i",
                    action = "store",
                    dest = "biominputfile",
                    metavar = "biomfile",
                    help = "[REQUIRED] sequence file",
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "outfilename",
                    metavar = "filename",
                    help = "[REQUIRED] outfile name",
                    required = True)
options = parser.parse_args()

#############################
# Import json formatted OTU #
#############################

import json
jsondata = open(options.biominputfile)
biom = json.load(jsondata)
jsondata.close()

from biom import Table
table = Table.from_json(biom)

print("")
print("Original OTU Table (without taxonomy)")
print("-------------------------------------")
print("")
print(table)
print("")

min_samplesize = int(min(table.sum(axis='sample')))
print("Subsampling to the smallest sample size: " + str(min_samplesize))

# Subsample
table_ss = table.subsample(min_samplesize)

# Output
outfile = open(options.outfilename, "w")
outfile.write("#OTU ID\t" + "\t".join(list(table_ss.ids("sample"))) + "\ttaxonomy" + "\n")
for otu in table_ss.ids("observation"):

    otu_counts = map(str, list(table_ss.data(otu, axis="observation")))
    taxonomy = table_ss.metadata(otu, "observation")["taxonomy"]
    confidence = table_ss.metadata(otu, "observation")["confidence"]
    
    outfile.write(otu + "\t" + "\t".join(otu_counts) + "\t" + "; ".join(taxonomy) + "\n")
        
exit(0)
