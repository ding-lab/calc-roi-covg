calcRoiCovg
===========
Given two BAM files, a reference sequence, and regions of interest, count the AT/CG/CpG sites with
sufficient read-depth in both BAMs.

Usage
-----
    calcRoiCovg <bam1> <bam2> <roi_file> <ref_seq_fasta> <output_file> [min_depth_bam1 min_depth_bam2 min_mapq]
    
    Defaults: min_depth_bam1 = 6, min_depth_bam2 = 8, min_mapq = 20
    NOTE: ROI file *must* be sorted by chromosome/contig names

This tool was originally designed to count base-pairs that have sufficient read-depth for variant
calling across two BAM files (case vs control). The base-pairs are further classified into AT, CG
(non-CpG), and CpG sites, with respect to the provided reference sequence. The resulting coverage
stats are reported for each region of interest (ROI).

Sample ROI files can be found under the 'data' subdirectory. Note that they use 1-based loci, and
*must* be sorted by chromosome or contig names.

Dependencies
------------
The Makefile assumes that you have the samtools source code in an environment variable $SAMDIR. If
you're working with a fresh install of Ubuntu, then simply follow these steps from your home dir:

    sudo apt-get install git libbam-dev zlib1g-dev
    mkdir src
    cd src
    git clone https://github.com/samtools/samtools.git
    export SAMDIR=$PWD/samtools
    git clone https://github.com/ckandoth/calc-roi-covg.git
    cd calc-roi-covg
    make
