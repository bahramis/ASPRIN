# ASPRIN
Allele-Specific Protein-RNA Interaction

Requirements
------------

ASPRIN is intended to be used in a Unix-based environment. It has been tested
on Mac OS and Linux.

Working installations of Python, pysam and samtools are required.

Python packages fisher, progressbar and mne are also required.

Installation
------------
All of the below installation instructions assume you have write access to the
installation directory for ASPRIN; if you want to install into a central
location for all users, you may need to prefix these with sudo.

### Configure ###
 Now type

    ./configure

This will configure the installation. ASPRIN is written partially in Python.
By default, it will use whichever python interpreter is invoked when you type
'Python' on the command line. If you wish to specify a different version of
Python, you can do this when running the ``configure`` script. Additionally,
if any of the required external software packages listed above are not on your
path, you can specify their location when running ``configure``. For example,
to specify the path to python and all of the exeternal software packages,
you would type the following:

    ./configure --with-python=/some/path/to/python

You needn't include them all, just the ones that cannot be found in your path.
The configure script will tell you if any requirements cannot be found, in
which case ASPRIN cannot be installed.

### Install ###
After configuration, ASPRIN is installed by typing:

    make install

This will install Python modules for whichever version of Python you are using.

Usage:
---------------------------------------

    $ asprin [-h] [-g Genotype file] [-c CLIP-seq mapped reads file]
              [-r RNA-seq mapped reads file] [-s dbSNP VCF file]
              [-e RADAR RNA editing database file]
              [-t Number of threads default: 25] [-a]

     ASPRIN: Allele Specific Protein-RNA Interaction

        optional arguments:
      -h, --help            show this help message and exit
      -s dbSNP VCF file
      -e RADAR RNA editing database file
      -t Number of threads (default: 20)
      -a                    Use all the variants

    required arguments:
      -g Genotype file
      -c CLIP-seq mapped reads file
      -r RNA-seq mapped reads file

To run ASPRIN, 3 arguments have to be provided"

1) The set of SNPs, in .vcf file format. This set of SNPs can either be 
genotype data obtained previously from any assay, or can be set of 
variants that are called from the RNA-seq data. This file is a vcf file

2) CLIP-seq data, in .bam file format, sorted based on coordinates and indexed.

3) RNA-seq data, in .bam file format, sorted based on coordinates and indexed.

Examples:
---------------------------------------

ASPRIN has 5 different modes of operation:


1) ASPRIN running on variants present in dbSNP only (SNPs)

```bash
    $ asprin -g example/HepG2_test_variants.vcf -s example/dbsnp_test.vcf
             -c example/RBFOX2_HepG2_clip_test.bam
             -r example/HepG2_total_rnaseq_test.bam
             -o output_snps.txt
```
2) ASPRIN running on variants present in RADAR only (RNA editing events)

```bash
    $ asprin -g example/HepG2_test_variants.vcf -e example/radar_test.txt
             -c example/RBFOX2_HepG2_clip_test.bam
             -r example/HepG2_total_rnaseq_test.bam
             -o output_editting.txt
```

3) ASPRIN running on variants in both dbSNP and RADAR 

```bash
    $ asprin -g example/HepG2_test_variants.vcf
             -s example/dbsnp_test.vcf -e example/radar_test.txt
             -c example/RBFOX2_HepG2_clip_test.bam
             -r example/HepG2_total_rnaseq_test.bam
             -o output_snps_and_editting.txt
```

4) ASPRIN running on all variants

```bash
    $ asprin -g example/HepG2_test_variants.vcf
             -c example/RBFOX2_HepG2_clip_test.bam
             -r example/HepG2_total_rnaseq_test.bam
             -o output_all.txt
```

5) ASPRIN running on all variants but label the SNPs and RNA editing events:

```bash
    $ asprin -g example/HepG2_test_variants.vcf
             -s example/dbsnp_test.vcf -e example/radar_test.txt
             -a
             -c example/RBFOX2_HepG2_clip_test.bam
             -r example/HepG2_total_rnaseq_test.bam
             -o output_all_labeled.txt
```

The example input and outputs are found in the folder called example.


dbSNP and RADAR databases:
---------------------------------------
dbSNP and RADAR under hg19 annotations databases can be downloaded from here:

dbSNP: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz

RADAR: http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19_v2.txt

For other versions of the databases please visit dbSNP and RADAR websites.

Enjoy!

Contacts and bug reports
------------------------
Emad Bahrami-Samani
ebs@ucla.edu

Yi Xing
yxing@ucla.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; we will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be.

Copyright and License Information
---------------------------------
Copyright (C) 2017 University of California, Los Angeles (UCLA)
Emad Bahrami-Samani, Yi Xing

Authors: Emad Bahrami-Samani, Yi Xing

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.
