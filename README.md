calc-roi-covg
=============
Given two BAM files, a reference sequence, and regions of interest, count the AT/CG/CpG sites with
sufficient read-depth in both BAMs.

Usage
-----
./calcRoiCovg <bam1> <bam2> <roi_file> <ref_seq_fasta> <output_file> [min_depth_bam1 min_depth_bam2 min_mapq]
Defaults: min_depth_bam1 = 6, min_depth_bam2 = 8, min_mapq = 20
NOTE: ROI file *must* be sorted by chromosome/contig names

This tool was designed for use with two BAM files (e.g. tumor and normal BAMs in cancer genomics)
that have sufficient read-depth for somatic variant calling. You can also use this to generate
similar coverage stats across a single BAM by passing the same BAM as bam1 and bam2.

Dependencies
------------
The Makefile assumes that you have the samtools source code in an environment variable $SAMDIR.
You can get that source from github:

    git clone https://github.com/samtools/samtools.git

If you're working with a fresh install of Ubuntu, then simply follow these steps from your home dir:

    sudo apt-get install git libbam-dev zlib1g-dev
    mkdir src
    cd src
    git clone https://github.com/samtools/samtools.git
    export SAMDIR=$PWD/samtools
    git clone https://github.com/ckandoth/calc-roi-covg.git
    cd calc-roi-covg
    make
