===========================================
gen-tools version 1.0 - README
Copyright 2009 Gavin Band, University of Oxford
===========================================

ACKNOWLEDGEMENTS:
================
This package includes the SNPHWE code from Wigginton et al, "A Note on Exact
Tests of Hardy-Weinberg Equilibrium", Am J Hum Genet (2005).

The package is built with the waf build system (included) 

DESCRIPTION:
============

This package comprises the qc-tool program for quality control of a set of SNP marker data,
together with some other utility programs.  The programs included are:

qc-tool: 		main quality control tool
gen-convert: 	convert between different gen file formats.
gen-compare: 	compare two gen files

There are also a number of tests (whose names begin with 'test') and some benchmark programs
(whose names beginning with 'benchmark').  

In addition, the following two programs are currently included but should really be in a seperate package:
generate-random-permutations-of-0-1-vector
gen-case-control-test


FILE FORMATS:
=============

For details of file formats used, see here:
http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html
and here:
http://www.well.ox.ac.uk/~gav/bgen_format.html

The programs identify gen file formats by the extension on the filename:
.gen
.gen.gz
.bgen
.bgen.gz

If none of the above, plain gen format is assumed.

PREREQUISITES
==============

gen-tools requires the following external libraries to be installed:

- Boost (version 1.36.1 or higher, I think.  I used 1.39.0 built using macports.)
- Zlib.

The configuration step (see below) will attempt to find these libraries.

COMPILATION
===========

Note: all built programs wind up the 'build/release' (release build) and 'build/default' (debug build) directories.

Compilation involves two steps: configuration and build.
To configure the package, from the top-level directory run
> ./waf-1.5.8 configure
To build the package, run
> ./waf-1.5.8

These commands build two versions of the programs: debug and release.  They can be found in the
build/default/ and build/release/ subdirectories respectively.

USAGE: qc-tool
==============

Given a collection of GEN files and a sample file, qc-tool can be used to:
- Print sample-wise statistics (missing data rate and heterozygosity) about the data (using -sample-stats flag).
- Filter out samples based on those statistics (using -sample-missing-rate and -heterozygosity).
- Print SNP-wise statistics (missing data rate, hwe exact test, minor allele frequency) about the data (using -snp-stats flag).
- Filter out SNPs based on those statistics (using -snp-missing-rate, -hwe, -maf)

You can also filter samples or snps based on inclusion or exclusion lists, a la gtool's select mode.

Other options:

-force :	Proceed even if there are warnings.

OPTIONS:

Data file options:
   -g <a>: Path of gen file to input.  Repeat this option, or use the numerical wildcard character '#', to specify several files.
  -og <a>: Path of gen file to output
  -os <a>: Path of sample file to output
   -s <a>: Path of sample file to input

Note that the -g option can take a filename containing a single '#' wildcard character.
This character matches any number from 1 to 100.  For example, if you have 23 files representing
data from 23 human chromosomes, these might be number file1.gen, file2.gen, ..., file23.gen you could run

> qc-tool -g file#.gen

Statistic file options:
  -sample-statistics <a>: Comma-seperated list of statistics to calculate in samplestat file.
       -sample-stats <a>: Output sample-wise statistics to the given file.
     -snp-statistics <a>: Comma-seperated list of statistics to calculate in genstat file.
          -snp-stats <a>: Output snp-wise statistics to the given file.

SNP filtering options:
               -hwe <a>: Filter out SNPs with HWE p-value less than or equal to the value specified.
           -maf <a> <b>: Filter out SNPs whose minor allele frequency lies outside the interval [a,b].
     -snp-excl-list <a>: Filter out SNPs whose SNP ID or RSID lies in the given file.
     -snp-incl-list <a>: Filter out SNPs whose SNP ID or RSID does not lie in the given file.
  -snp-interval <a> <b>: Filter out SNPs with position outside the interval [a,b].
  -snp-missing-rate <a>: Filter out SNPs with missing data rate greater than or equal to the value specified.

Sample filtering options:
   -heterozygosity <a> <b>: Filter out samples with heterozygosity outside the inteval [a,b].
     -sample-excl-list <a>: Filter out samples whose sample ID lies in the given file.
     -sample-incl-list <a>: Filter out samples whose sample ID does not lie in the given file.
  -sample-missing-rate <a>: Filter out samples with missing data rate greater than the value specified.

Other options:
  -diagnose-sample-filter: Print diagnostic information about each filtered out sample.
     -diagnose-snp-filter: Print diagnostic information about each filtered out snp.
                   -force: Ignore warnings and proceed with requested action.

Note: qc-tool will emit a warning and quit if it thinks the options supplied don't make sense.  To override this, use the -force option.


USAGE: gen-convert
==================

gen-convert converts one or more input gen files into one or more output gen files (possibly in a different gen file format).
Typical usage is:

> gen-convert -g (input files) -og (output files)

The -g option must appear the same number of times as the -og option.
As with qc-tool, a '#' wildcard character may appear in the input file and matches a number from 1 to 100.
If such a wildcard appears in the -g option, a similar wildcard may appear for the -og option.  For example,
with the files as above we could write

> gen-convert -g file#.gen -og file#.bgen

The result would be 23 files names file1.bgen, file2.bgen, ... each one a bgen format version of the 
corresponding input file.


USAGE: gen-compare
==================

Compare two gen (or bgen, gen.gz, etc. files).  This is a simple program which will print out mismatches between 
two files.  It can be used to verify the output of gen-convert.  For example:

> gen-compare -g1 file#.gen -g2 file#.bgen


Thanks,
Gavin Band.
