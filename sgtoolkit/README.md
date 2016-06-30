SGToolKit
====

Command line utilities for illumina seqeunce data processing. These tools are for our internal use.

Installing from source: 

```sh
wget https://github.com/hsgweon/sgtoolkit/archive/master.zip
cd sgtoolkit-master
python setup.py clean --all
python setup.py install --prefix=$HOME/sgtoolkit
```

Make sure you have set `$PATH` to point to `$HOME/sgtoolkit`

SGTK is an on-going project and written for internal projects. For usage information, type commands with `-h`

SGTK Pipeline consists of three steps.

```sh
sgtk_make_file_pairs_list.py -i rawdata_dir
sgtk_prepseqs.py -i rawdata_dir -o sgtk_prepseqs -l file_pairs_list.txt
sgtk_processseqs.py -i sgtk_prepseqs/prepped.fasta -r 16S -l file_pairs_list.txt -o sgtk_processseq
```

sgtk_make_file_pairs_list.py
----------------------------

This will try to deduce sample_ids and their corresponding forward read filename, reverse read filename from a directory with raw reads from Illumina sequencing runs. It was first intended to be used for PIPITS.

Example:

```sh
sgtk_make_file_pairs_list.py -i rawdata_dir
```

For help, type

```sh
sgtk_make_file_pairs_list.py -h
```

sgtk_prepseqs.py
-----------

By default, Illumina sequencers generate demultiplexed sequence files. So all you need is a directory with all those paired fastq files to run SGTK which will reindex, join, quality filter, convert and merge your demultiplexed reads. The resulting file can then be processed by PIPITS, QIIME etc.

*Dependencies*

-   PEAR (<http://sco.h-its.org/exelixis/web/software/pear>) - N.B. PEAR
    prohibits commercial use of the code. See its page for detail.
-   FASTX-Toolkit (<http://hannonlab.cshl.edu/fastx_toolkit>)
    *Available as a Bio-Linux package*

Example:

```sh 
sgtk_prepseqs.py -i rawdata_dir -o sgtk_prepseqs -l file_pairs_list.txt
```

For help, type

```sh
sgtk_prepseqs.py -h
```

sgtk_processseqs.py
-----------

Pipeline for 16S & 18S amplicon analysis. Intended for our internal use in Centre for Ecology and Hydrology. Not rigorously tested, and is continually being updated.

*Dependencies*

-   BIOM-FORMAT v.2.x
    (<https://pypi.python.org/pypi/biom-format/1.3.1>)
   
    *Available as a Bio-Linux package*

-   VSEARCH (<https://github.com/torognes/vsearch>)
   
    *Available as a Bio-Linux package*

-   RDP Classifier 2.9 or above
    (<http://sourceforge.net/projects/rdp-classifier>) - N.B. RDP
    Classifier comes with a jar file.

-   QIIME (<http://qiime.org>)
   
    *Available as a Bio-Linux package*

Example:

```sh
sgtk_processseqs.py -i sgtk_prepseqs/prepped.fasta -r 16S -l file_pairs_list.txt -o sgtk_processseq
```

For detailed information on options:

```sh
sgtk_processseqs.py -h
```
