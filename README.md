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

Install
-------
The Makefile assumes that you have the samtools source code in an environment variable `$SAMDIR`. If
you don't know what that means, then simply follow these steps from any directory that you have
permissions to write into:

Install some prerequisite packages if you are using Debian or Ubuntu:

    sudo apt-get install git libbam-dev zlib1g-dev

If you are using Fedora, CentOS or RHEL, you'll need these packages instead:

    sudo yum install git samtools-devel zlib-devel

Clone the samtools and calc-roi-covg repos, and build the `calcRoiCovg` binary:

    git clone https://github.com/samtools/samtools.git
    export SAMDIR=$PWD/samtools
    git clone https://github.com/ding-lab/calc-roi-covg.git
    cd calc-roi-covg
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have su permissions, then
I recommend dumping it in the system directory for locally compiled packages:

    sudo mv calcRoiCovg /usr/local/bin/
