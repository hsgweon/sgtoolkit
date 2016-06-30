#!/usr/bin/python

"""
PROCESSSEQ = Generate an OTU Table from PREPSEQ'ed sequences

This works with demultiplexed files from Illumina platform.
"""

import argparse, sys, os, argparse, shutil, subprocess, textwrap, logging, gzip, bz2

__author__ = "Hyun Soon Gweon"
__copyright__ = "Hyun Soon Gweon"
__credits__ = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__ = "GPL"
__maintainer__ = "Hyun Soon Gweon"
__email__ = "hyugwe@ceh.ac.uk"

Red = '\033[91m'
Green = '\033[92m'
Blue = '\033[94m'
Cyan = '\033[96m'
White = '\033[97m'
Yellow = '\033[93m'
Magenta = '\033[95m'
Grey = '\033[90m'
Black = '\033[90m'
Default = '\033[0m'

QIIME = "qiime"
BIOM = "biom"
VSEARCH = "vsearch"

#DB_16S_CHIMERA = "/home/hyugwe/shared/db/greengenes/gg_13_5/gg_13_5_otus/rep_set/97_otus.fasta"
#DB_16S_ASSIGNMENT_REF_FASTA = "/home/hyugwe/shared/db/greengenes/gg_13_5/gg_13_5.fasta"
#DB_16S_ASSIGNMENT_REF_TAXONOMY = "/home/hyugwe/shared/db/greengenes/gg_13_5/gg_13_5_taxonomy.txt"
#DB_16S_ALIGNMENT = "/home/hyugwe/shared/db/greengenes/gg_13_5/gg_13_5_otus/rep_set_aligned/85_otus.fasta"

#DB_18S_CHIMERA = "/home/hyugwe/shared/db/silva/SILVA_119/Silva119_release/rep_set_eukaryotes/99/Silva_119_rep_set99_18S.fna"
#DB_18S_ASSIGNMENT_REF_FASTA = "/home/hyugwe/shared/db/silva/SILVA_119/Silva119_release/rep_set_eukaryotes/99/Silva_119_rep_set99_18S.fna"
#DB_18S_ASSIGNMENT_REF_TAXONOMY = "/home/hyugwe/shared/db/silva/SILVA_119/Silva119_release/taxonomy_eukaryotes/99/taxonomy_99_7_levels_18S.txt"
#DB_18S_ALIGNMENT = "/home/hyugwe/shared/db/silva/SILVA_119/Silva119_release/core_alignment/core_Silva119_alignment.fna"

DB_18SPR2_CHIMERA = "/home/hyugwe/shared/db/silva/SILVA_119/Silva119_release/rep_set_eukaryotes/99/Silva_119_rep_set99_18S.fna"
DB_18SPR2_ASSIGNMENT_REF_FASTA = "/home/hyugwe/shared/db/prr/mothur_qiime_gb203.fasta"
DB_18SPR2_ASSIGNMENT_REF_TAXONOMY = "/home/hyugwe/shared/db/prr/qiime_gb203_taxo.txt"
DB_18SPR2_ALIGNMENT = "/home/hyugwe/shared/db/silva/SILVA_119/Silva119_release/core_alignment/core_Silva119_alignment.fna"


VSEARCH_CLUSTER_THRESHOLD = "0.97"

DEREP_DIR = ""
CLUSTER_DIR = ""
OTUWITHTAXONOMY_DIR = ""
REMOVECHIMERA_DIR = ""

REMOVEUNMATCHEDSEQUENCES_DIR = ""

REMAPPED_DIR = ""
UC2OTUTABLE_DIR = ""
CLASSIFYREPSET_DIR = ""


def run_cmd(command):
    
    FNULL = open(os.devnull, 'w')
    if options.verbose:
        print(Green + command + Default)
        p = subprocess.Popen(command, shell=True)
    else:
        p = subprocess.Popen(command, shell=True, stdout=FNULL, stderr=FNULL)
    p.wait()
    FNULL.close()
    if p.returncode != 0:
        print("None zero returncode: " + command)
        exit(1)

def getsamplelistfromfasta(options):

    global DEREP_DIR

    cmd = " ".join(["sgtk_getsamplelistfromfasta.py",
                    "-i", options.inputfasta,
                    "-o", options.outputdir + "/samplelist.lst"])

    print(Blue + "Generating a sample list from the input sequences..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)
        if os.stat(options.outputdir + "/samplelist.lst").st_size == 0:
            print("Sample list file is empty. Processing stopped.")
            exit(0)

def derep(options):
    
    global DEREP_DIR

    # Derep with sgtk
    if not options.printonly:
        if os.path.exists(DEREP_DIR):
            shutil.rmtree(DEREP_DIR)
        os.mkdir(DEREP_DIR)

    minuniquesize = 2 # 2 is the default. I.e. remove unique sequences.
    if options.includeuniqueseqs:
        minuniquesize = 1
    cmd = " ".join([VSEARCH, "--derep_fulllength", options.inputfasta, 
                    "--output", DEREP_DIR + "/input_dereplicated.fasta", 
                    "--minuniquesize", str(minuniquesize), 
                    "--sizeout",
                    "--threads", options.threads])


    print(Blue + "Dereplicating and removing unique sequences..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)

    # Check if the file is empty
        if os.stat(DEREP_DIR + "/input_dereplicated.fasta").st_size == 0:
            print("After dereplicating and removing unique sequences, there aren't zero sequences! Processing stopped.")
            exit(0)


def cluster(options):

    global DEREP_DIR
    global CLUSTER_DIR

    if not options.printonly:
        if os.path.exists(CLUSTER_DIR):
            shutil.rmtree(CLUSTER_DIR)
        os.mkdir(CLUSTER_DIR)

    cmd = " ".join([VSEARCH,
                    "--cluster_fast", DEREP_DIR + "/input_dereplicated.fasta",
                    "--id", VSEARCH_CLUSTER_THRESHOLD,
                    "--uc", CLUSTER_DIR + "/all_clusters.uc",
                    "--centroids", CLUSTER_DIR + "/centroids.fasta",
                    "--thread", options.threads])

    print(Blue + "Clustering sequences into OTUs..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


def removeChimera(options):

    global REMOVECHIMERA_DIR
    global CLUSTER_DIR

    if not options.printonly:
        if os.path.exists(REMOVECHIMERA_DIR):
            shutil.rmtree(REMOVECHIMERA_DIR)
        os.mkdir(REMOVECHIMERA_DIR)

    '''
    if options.region == "16S":
        DB_CHIMERA = DB_16S_CHIMERA
    elif options.region == "18S":
        DB_CHIMERA = DB_18S_CHIMERA

    cmd = " ".join([VSEARCH,
                    "--uchime_ref", CLUSTER_DIR + "/centroids.fasta",
                    "--db", DB_CHIMERA,
                    "--notrunclabels",
                    "--nonchimeras", REMOVECHIMERA_DIR + "/centroids_chimeraless.fasta",
                    "--chimeras", REMOVECHIMERA_DIR + "/centroids_chimeras.fasta"])
    '''

    # de novo
    cmd = " ".join([VSEARCH,
                    "--uchime_denovo", CLUSTER_DIR + "/centroids.fasta",
                    "--notrunclabels",
                    "--nonchimeras", REMOVECHIMERA_DIR + "/centroids_chimeraless.fasta",
                    "--chimeras", REMOVECHIMERA_DIR + "/centroids_chimeras.fasta"])

    print(Blue + "Detecting and removing chimerass..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


def removeUnmatchedSequences(options):

    global REMOVEUNMATCHEDSEQUENCES_DIR

    if not options.printonly:
        if os.path.exists(REMOVEUNMATCHEDSEQUENCES_DIR):
            shutil.rmtree(REMOVEUNMATCHEDSEQUENCES_DIR)
        os.mkdir(REMOVEUNMATCHEDSEQUENCES_DIR)

    if options.region == "16S":
        DB_REF = "$DB_16S_ASSIGNMENT_REF_FASTA"
    elif options.region == "18S":
        DB_REF = "$DB_18S_ASSIGNMENT_REF_FASTA"

    cmd = " ".join([VSEARCH,
                    "--usearch_global", REMOVECHIMERA_DIR + "/centroids_chimeraless.fasta",
                    "--db", DB_REF,
                    "--id 0.70",
                    "--matched", REMOVEUNMATCHEDSEQUENCES_DIR +    "/centroids_chimeraless_matched.fasta",
                    "--notmatched", REMOVEUNMATCHEDSEQUENCES_DIR + "/centroids_chimeraless_notmatched.fasta",
                    "--threads", options.threads])

    print(Blue + "Removing sequences below a similarity threshold (70%) against known sequences (Greengenes for 16S, Silva for 18S)..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)



def renameRepset(options):

    global REMOVEUNMATCHEDSEQUENCES_DIR

    print(Blue + "Re-indexing Representative Sequences (OTUs)..." + Default)
    if not options.printonly:
        
        handle_in = open(REMOVEUNMATCHEDSEQUENCES_DIR + "/centroids_chimeraless_matched.fasta", "rU")
        handle_out = open(REMOVEUNMATCHEDSEQUENCES_DIR + "/centroids_chimeraless_matched_reindexed.fasta", "w")
        
        counter = 1
        
        for line in handle_in:
            if line.startswith(">"):
                newlabel = line[1:].split(";")[0]
                handle_out.write(">OTU" + str(counter) + "\n")
                counter += 1
            else:
                handle_out.write(line.rstrip() + "\n")
        handle_in.close()
        handle_out.close()

    # Copy repseqs to main dir
    cmd = " ".join(["cp -r", REMOVEUNMATCHEDSEQUENCES_DIR + "/centroids_chimeraless_matched_reindexed.fasta", options.outputdir + "/repseqs.fasta"])
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


def mapReadsOntoRepset(options):

    global REMAPPED_DIR

    if not options.printonly:
        if os.path.exists(REMAPPED_DIR):
            shutil.rmtree(REMAPPED_DIR)
        os.mkdir(REMAPPED_DIR)

    cmd = " ".join([VSEARCH, 
                    "--usearch_global", options.inputfasta, 
                    "--db", REMOVEUNMATCHEDSEQUENCES_DIR + "/centroids_chimeraless_matched_reindexed.fasta", 
                    "--id", VSEARCH_CLUSTER_THRESHOLD, 
                    "--uc", REMAPPED_DIR + "/otus.uc",
                    "--threads", options.threads])

    print(Blue + "Map reads onto OTUs..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


def uc2otutable(options):
    
    global UC2OTUTABLE_DIR

    if not options.printonly:
        if os.path.exists(UC2OTUTABLE_DIR):
            shutil.rmtree(UC2OTUTABLE_DIR)
        os.mkdir(UC2OTUTABLE_DIR)

    cmd = " ".join(["sgtk_uc2otutable.py",
                    "-i", REMAPPED_DIR + "/otus.uc",
                    "-o", UC2OTUTABLE_DIR + "/otu_table.txt",
                    "-l", options.outputdir + "/samplelist.lst"])

    print(Blue + "Converting uc file into an OTU table..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


    # Convert to biom
    try:
        os.remove(UC2OTUTABLE_DIR + "/otu_table.biom")
    except OSError:
        pass
    cmd = " ".join([BIOM,
                    "convert", 
                    "-i", UC2OTUTABLE_DIR + "/otu_table.txt", 
                    "-o", UC2OTUTABLE_DIR + "/otu_table.biom", 
                    "--table-type=\"OTU table\" --to-json"])

    print(Blue + "Converting OTU format..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


def classifyRepset(options):

    global CLASSIFYREPSET_DIR

    if not options.printonly:
        if os.path.exists(CLASSIFYREPSET_DIR):
            shutil.rmtree(CLASSIFYREPSET_DIR)
        os.mkdir(CLASSIFYREPSET_DIR)

    if options.region == "16S":
        DB_ASSIGNMENT_REF_FASTA = "$DB_16S_ASSIGNMENT_REF_FASTA"
        DB_ASSIGNMENT_REF_TAXONOMY = "$DB_16S_ASSIGNMENT_REF_TAXONOMY"
    elif options.region== "18S":
        DB_ASSIGNMENT_REF_FASTA = "$DB_18S_ASSIGNMENT_REF_FASTA"
        DB_ASSIGNMENT_REF_TAXONOMY = "$DB_18S_ASSIGNMENT_REF_TAXONOMY"

    # RDP Classifier
    cmd = " ".join(["qiime", "assign_taxonomy.py", 
                    "-i", REMOVEUNMATCHEDSEQUENCES_DIR + "/centroids_chimeraless_matched_reindexed.fasta", 
                    "-r", DB_ASSIGNMENT_REF_FASTA,
                    "-t", DB_ASSIGNMENT_REF_TAXONOMY,
                    "-m", "rdp",
                    "--rdp_max_memory 15000",
                    "-o", CLASSIFYREPSET_DIR])

    print(Blue + "Assigning taxonomy with RDP Classifier..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


    '''

    # BLAST
    cmd = " ".join(["qiime", "assign_taxonomy.py",
                    "-i", REMOVECHIMERA_DIR + "/centroids_chimeraless_reindexed.fasta",
                    "-r", DB_ASSIGNMENT_REF_FASTA,
                    "-t", DB_ASSIGNMENT_REF_TAXONOMY,
                    "-m", "blast",
                    "-o", CLASSIFYREPSET_DIR + "_blast"])

    print(Blue + "Assigning taxonomy with BLAST..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)
    '''

def addTaxonomyToOTU(_input_dir, _input_otu_biom, _output_dir, _output_otu_biom, _output_otu_tabular, _tax_assignments_dir, _tax_assignments, _printonly):

    if not _printonly:
        if os.path.exists(_output_dir):
            shutil.rmtree(_output_dir)
        os.mkdir(_output_dir)
    
    # Adding RDP_CLASSIFIER output to OTU table
    cmd = " ".join([BIOM, "add-metadata", 
                    "-i", _input_dir + "/" + _input_otu_biom, 
                    "-o", _output_dir + "/" + _output_otu_biom, 
                    "--observation-metadata-fp", _tax_assignments_dir + "/" + _tax_assignments,
                    "--observation-header", "OTUID,taxonomy,confidence", 
                    "--sc-separated", "taxonomy", 
                    "--float-fields", "confidence",
                    "--output-as-json"])
    print(Blue + "Adding assignment to OTU Table..." + Default)
    if _printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)

    cmd = " ".join([BIOM, "convert",
                    "-i", _output_dir + "/" + _output_otu_biom,
                    "-o", _output_dir + "/" + _output_otu_tabular,
                    "--header-key taxonomy",
                    "--to-tsv"])
    print(Blue + "Converting BIOM to TXT..." + Default)
    if _printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)

    '''
    cmd = " ".join(["cp -r", _output_dir + "/" + _output_otu_biom, ". ;",
                    "cp -r", _output_dir + "/" + _output_otu_tabular, ". ;"])
    if _printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)
    '''

def summarizeTable(options):

    # Summarize Table
    cmd = " ".join([BIOM, "summarize-table",
                    "-i", OTUWITHTAXONOMY_DIR + "/otu_table.biom",
                    "-o", OTUWITHTAXONOMY_DIR + "/summary_otu_table.txt"])

    print(Blue + "Summarising OTU table..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)



def summarizeTaxa(options):

    cmd = " ".join([QIIME, "summarize_taxa.py",
                    "-i", OTUWITHTAXONOMY_DIR + "/otu_table.biom",
                    "-o", OTUWITHTAXONOMY_DIR + "/taxa_summary",
                    "-L 2,3,4,5,6,7 -a"])

    print(Blue + "Summarising taxonomic information..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)
    

    plot_taxa_summary_input = []
    for i in range(2, 8):
        plot_taxa_summary_input.append(OTUWITHTAXONOMY_DIR + "/taxa_summary/otu_table_L" + str(i) + ".txt")
    cmd = " ".join([QIIME, "plot_taxa_summary.py",
                    "-i", ",".join(plot_taxa_summary_input),
                    "-o", OTUWITHTAXONOMY_DIR + "/taxa_summary/taxa_summary_plots"])

    print(Blue + "Plotting some basic plots..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)

    # Copy files
    cmd = " ".join(["cp -r", OTUWITHTAXONOMY_DIR + "/taxa_summary", options.outputdir + ";",
                    "cp -r", OTUWITHTAXONOMY_DIR + "/summary_otu_table.txt", options.outputdir + ";"
                    ])
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


def build_phylogenetic_tree(options):

    global PHYLOGENY_DIR

    if not options.printonly:
        if os.path.exists(PHYLOGENY_DIR):
            shutil.rmtree(PHYLOGENY_DIR)
        os.mkdir(PHYLOGENY_DIR)

    if options.region == "16S":
        DB_ALIGNMENT = "$DB_16S_ALIGNMENT"
    elif options.region== "18S":
        DB_ALIGNMENT = "$DB_18S_ALIGNMENT"

    # Align
    cmd = " ".join([QIIME, "align_seqs.py",
                    "-i", REMOVEUNMATCHEDSEQUENCES_DIR + "/centroids_chimeraless_matched_reindexed.fasta",
                    "-o", PHYLOGENY_DIR,
                    "-t", DB_ALIGNMENT,
                    "-p", "0.01",
                    "-e", "0"])

    print(Blue + "Aligning representative sequences..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)
        
    # Build a phylogetic tree
    cmd = " ".join([QIIME, "make_phylogeny.py",
                    "-i", PHYLOGENY_DIR + "/centroids_chimeraless_matched_reindexed_aligned.fasta",
                    "-o", PHYLOGENY_DIR + "/repseqs.tre"])
    
    print(Blue + "Building a phylogenetic tree..." + Default)
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)

    # Copy files
    cmd = " ".join(["cp -r", PHYLOGENY_DIR + "/repseqs.tre", options.outputdir])
    if options.printonly:
        print(Green + cmd + Default)
    else:
        run_cmd(cmd)


def removeIntermediateFiles(options):

    if options.printonly:
        pass
    else:
        if not options.keep:
            print(Blue + "Removing intermediate files..." + Default)
            shutil.rmtree(DEREP_DIR)
            shutil.rmtree(CLUSTER_DIR)
            shutil.rmtree(REMOVECHIMERA_DIR)
            shutil.rmtree(REMOVEUNMATCHEDSEQUENCES_DIR)
            shutil.rmtree(REMAPPED_DIR)
            shutil.rmtree(UC2OTUTABLE_DIR)
            shutil.rmtree(CLASSIFYREPSET_DIR)
            shutil.rmtree(OTUWITHTAXONOMY_DIR)
            shutil.rmtree(PHYLOGENY_DIR)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Processes PREPSEQ'ed sequences and generates an OTU table. For CEH internal usage only for the moment. Contact Soon Gweon (hyugwe@ceh.ac.uk) for any enquiry. Currently supports 16S and 18S. Use PIPITS for ITS.")
    parser.add_argument(
        "-i",
        action = "store",
        dest = "inputfasta",
        metavar = "<FASTA>",
        help = "[REQUIRED] Prepped sequences (preferably by PREPSEQ)",
        required = True)
    parser.add_argument(
        "--region",
        action = "store",
        dest = "region",
        help = "[REQUIRED] region: \"16S\", \"18S\"",
        required = True,
        choices = ["16S", "18S"])
    parser.add_argument(
        "-o",
        action = "store",
        dest = "outputdir",
        metavar = "<DIR>",
        help = "Directory to output results [default: sgtk_processseqs]",
        default = "sgtk_processseqs",
        required = False)
    parser.add_argument(
        "--includeuniqueseqs",
        action = "store_true",
        dest = "includeuniqueseqs",
        help = "[REQUIRED] PIPITS by default removes unique sequences before clustering. This means you wouldn't have any singletons. If you want singletons, then choose this option. It can take much longer to process.",
        required = False)
    parser.add_argument(
        "-t",
        action = "store",
        dest = "threads",
        metavar = "<INT>",
        help = "Number of Threads [default: 1]",
        default = "1",
        required = False)
    parser.add_argument(
        "-p, --printonly",
        action = "store_true",
        dest = "printonly",
        help = "Do NOT execute commands, but just print commands",
        required = False)
    parser.add_argument(
        "-v, --verbose",
        action = "store_true",
        dest = "verbose",
        help = "Verbose mode",
        required = False)
    parser.add_argument(
        "-r, --retain",
        action = "store_true",
        dest = "keep",
        help = "Keep intermediate files. These are laaaaarge, so unless there is a good reason, please don't choose this option.",
        required = False)
    options = parser.parse_args()

    DEREP_DIR = options.outputdir + "/1_dereplicated"
    CLUSTER_DIR = options.outputdir + "/2_repset_preliminary"
    REMOVECHIMERA_DIR = options.outputdir + "/3_repset"
    REMOVEUNMATCHEDSEQUENCES_DIR = options.outputdir + "/4_remove_unmatched_sequences"
    REMAPPED_DIR = options.outputdir + "/5_remapped"
    UC2OTUTABLE_DIR = options.outputdir + "/6_otu_table_preliminary"
    CLASSIFYREPSET_DIR = options.outputdir + "/7_taxonomic_assignment"
    OTUWITHTAXONOMY_DIR = options.outputdir + "/8_otu_table"
    PHYLOGENY_DIR = options.outputdir + "/9_phylogenetic_tree"


    if not options.printonly:
        if os.path.exists(options.outputdir):
            shutil.rmtree(options.outputdir)
            os.mkdir(options.outputdir)
            pass
        else:
            os.mkdir(options.outputdir)


    getsamplelistfromfasta(options)
    derep(options)
    cluster(options)
    removeChimera(options)
    removeUnmatchedSequences(options)
    renameRepset(options)
    mapReadsOntoRepset(options)
    uc2otutable(options)
    classifyRepset(options)


    # RDP
    addTaxonomyToOTU(_input_dir = UC2OTUTABLE_DIR, 
                     _input_otu_biom = "otu_table.biom", 
                     _output_dir = OTUWITHTAXONOMY_DIR, 
                     _output_otu_biom = "otu_table.biom", 
                     _output_otu_tabular = "otu_table.txt", 
                     _tax_assignments_dir = CLASSIFYREPSET_DIR, 
                     _tax_assignments = "centroids_chimeraless_matched_reindexed_tax_assignments.txt", 
                     _printonly = options.printonly)

    '''
    # BLAST
    addTaxonomyToOTU(_input_dir = UC2OTUTABLE_DIR,
                     _input_otu_biom = "otu_table.biom",
                     _output_dir = OTUWITHTAXONOMY_DIR + "_BLAST",
                     _output_otu_biom = "otu_table.biom",
                     _output_otu_tabular = "otu_table.txt",
                     _tax_assignments_dir = CLASSIFYREPSET_DIR + "_blast",
                     _tax_assignments = "centroids_chimeraless_reindexed_tax_assignments.txt",
                     _printonly = options.printonly)
    '''

    summarizeTable(options)
    build_phylogenetic_tree(options)
    summarizeTaxa(options)
    removeIntermediateFiles(options)

    exit(0)
