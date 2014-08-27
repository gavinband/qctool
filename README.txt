===========================================
qctool - README
Copyright 2009-2011 Gavin Band, University of Oxford
===========================================

ACKNOWLEDGEMENTS:
================
This package includes the SNPHWE code from Wigginton et al, "A Note on Exact
Tests of Hardy-Weinberg Equilibrium", Am J Hum Genet (2005).

The package is built with the waf build system (included with the source).  See http://code.google.com/p/waf/

This package includes a VCF format implementation supplied by Quang Le, found under genfile/test/CVCFT.

DESCRIPTION:
============

This package comprises the QCTOOL program for quality control of a set of SNP marker data,
together with some other utility programs.  The programs included are:

qctool: 		main quality control tool.
inthinnerator:  thin SNPs.
gen-grep:		find entries in a gen file (or files).

FILE FORMATS:
=============

QCTOOL supports the following file formats.
* GEN (or IMPUTE) format, possibly gzipped.  See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html.
* BGEN format. See http://www.well.ox.ac.uk/~gav/bgen_format.html.
* VCF format v4.0 and above.  See http://vcftools.sourceforge.net/specs.html.
* SHAPEIT and IMPUTE-style haplotype formats
* Custom formats

The programs identify file formats by inspecting the extension on the filename:
.gen
.gen.gz
.bgen
.vcf
.vcf.gz

If none of the above, plain gen format is assumed.  The options -[i|o]filetype can be used to specify the input/output file types.

PREREQUISITES
==============

QCTOOL requires the following external libraries:

- zlib (http://zlib.net/)

QCTOOL contains the source code to boost.

COMPILATION
===========

Compilation involves two steps: configuration and build.
To configure the package, from the top-level directory run
> ./waf-1.5.18 configure
To build the package, run
> ./waf-1.5.18

Then main qctool executable can be found under build/release/.  For debugging purposes, a debug build is also built under build/default/.

USAGE: qctool
==============

Given a collection of GEN files and a sample file, qctool can be used to:
- Print sample-wise statistics (missing data rate and heterozygosity) about the data (using -sample-stats flag).
- Filter out samples based on those statistics (using -sample-missing-rate and -heterozygosity).
- Print SNP-wise statistics (missing data rate, hwe exact test, minor allele frequency, information) about the data (using -snp-stats flag).
- Filter out SNPs based on those statistics (using -snp-missing-rate, -hwe, -maf, -info)
- Align SNPs to a reference genome using a strand file
- Match SNPs between cohorts, generating an overlap set of SNPs.
- Infer relatedness between individuals
- Manipulate covariates and phenotypes present in the sample file.

You can also filter samples or snps based on inclusion or exclusion lists, a la gtool's select mode.

Note: qctool will emit a warning and quit if it thinks the options supplied don't make sense.  To override this, use the -force option.

Here is a full list of options.

Usage: qctool <options>

OPTIONS:
Input file options:
                                -g <a>: Path of gen file(s) to input.  The given filename may contain the wildcard chara-
                                        cter '#', which expands to match a two-character chromosome identifier.  (For ex-
                                        ample, "qctool -g myfile_#.gen" will find all files of the form "myfile_01.gen",
                                        "myfile_02.gen", etc.)  Only Human autosomes are matched this way.
                                        This option may be repeated, in which case each invocation is treated as a seper-
                                        ate cohort and cohortsare joined together to create one big dataset.
                                -s <a>: Path of sample file to input.  If specified, this option must occur as often as 
                                        the -g option to specify one sample file per cohort.

Sample exclusion options:
                     -excl-samples <a>: Filter out samples whose sample ID lies in the given file.
                     -incl-samples <a>: Filter out samples whose sample ID does not lie in the given file.

Options for adjusting sample data:
                     -missing-code <a>: Specify a comma-separated list of strings to be treated as missing values when e-
                                        ncountered in the sample file(s).  Defaults to "NA".
               -quantile-normalise <a>: Quantile normalise each specified continuous phenotype or covariate by ranking i-
                                        ts values and mapping to quantiles of the standard normal distribution N(0,1). T-
                                        ies are handled by sending tied values to the average of the corresponding quant-
                                        iles.The argument should be a comma-separated list of column names from the samp-
                                        le file.

SNP exclusion options:
                       -excl-rsids <a>: Exclude all SNPs whose RSID is in the given file(s) from the analysis.
                      -excl-snpids <a>: Exclude all SNPs whose SNPID is in the given file(s) from the analysis.
               -excl-snps-matching <a>: Filter out snps whose rsid or SNPID matches the given value. The value should be
                                        a string which can contain a % wildcard character (which matches any substring).
                                        If you use *, you should place the argument in quotes.Optionally, prefix the arg-
                                        ument with snpid~ or rsid~ to only match against the SNPID or rsid fields.
                       -incl-rsids <a>: Exclude all SNPs whose RSID is not in the given file(s) from the analysis.
                      -incl-snpids <a>: Exclude all SNPs whose SNPID is not in the given file(s) from the analysis.
               -incl-snps-matching <a>: Filter out snps whose rsid or SNPID does not match the given value. The value sh-
                                        ould be a string which can contain a % wildcard character (which matches any sub-
                                        string). Optionally, prefix the argument with snpid~ or rsid~ to only match agai-
                                        nst the SNPID or rsid fields.

Options for adjusting SNPs:
                -assume-chromosome <a>: Treat each SNP whose chromosome cannot be determined as though it lies on the sp-
                                        ecified chromosome.
             -match-alleles-to-cohort1: Specify that alleles (and corresponding genotypes) in all cohorts should be swit-
                                        ched, if necessary, so as to match the alleles of the first cohort.  This does n-
                                        ot perform allele complementation, but you can use the -strand option to complem-
                                        ent alleles first.
                 -snp-match-fields <a>: By default, matching SNPs between cohorts uses all the available fields (positio-
                                        n, rsid, SNPID, and alleles.) Use this option to specify a comma-separated subse-
                                        t of those fields to use instead. The first entry must be "position". This optio-
                                        n can be used, for example, when cohorts are typed on different platforms so hav-
                                        e different SNPID fields.  Defaults to "position,rsid,SNPID,alleles".
                           -strand <a>: Path of strand file(s) to input.  If specified, this option must occur the same 
                                        number of times as the -g option, to specify one intensity file per cohort.
          -translate-snp-positions <a>: Specify a "dictionary" of chromosome / position to chromosome / position mapping-
                                        s. (This should come as a four-column file with source_chromosome source_positio-
                                        n target_chromosome and target_position columns.) Positions of SNPs will be mapp-
                                        ed through this dictionary before processing.

Output file options:
                               -og <a>: Override the auto-generated path(s) of the output gen file(s) to use when filter-
                                        ing.  (By default, the paths are formed by adding ".fltrd" to the input gen file-
                                        name(s).)  If this option is supplied, it must appear the same number of times a-
                                        s the -g option. If the corresponding occurence of -g uses a '#' wildcard charac-
                                        ter, the '#' character can also be used here to specify numbered output files co-
                                        rresponding to the input files.
                               -os <a>: Override the auto-generated path of the output sample file.  

Pedigree file options:
                               -ip <a>: Input a pedigree from the specified file. The first six columns of this file sho-
                                        uld represent a PED format pedigree, according to the spec on the PLINK website.
                                        (Other columns are ignored.) Ids are treated as non-whitespace strings and sex c-
                                        an be either "1" or "M" (male) or "2" or "F" (female) or "other".
                               -op <a>: Output a pedigree file instead of a GEN-type file. You must also input a pedigre-
                                        e using -ip for this to work.

VCF file options:
               -vcf-genotype-field <a>: Specify the name of the field in a VCF or VCF-like file from which genotype call
                                        probabilities will be taken. This must either be the GT field (in which case cal-
                                        ls are treated as certain calls with the given genotype), or a field of type Flo-
                                        at with Number equal to "3", "G", or the missing value ".".  Defaults to "GT".

Statistic calculation options:
                         -sample-stats: Calculate and output sample-wise statistics.
             -sample-stats-columns <a>: Comma-seperated list of statistics to output in the sample-wise statistics file.
                                        By default, the columns are: ID1, ID2, missing, and heterozygosity.  Defaults to
                                        "ID1, ID2, missing, heterozygosity".
                -sample-stats-file <a>: Override the auto-generated path of the file in which sample-wise statistics wil-
                                        l be output.
                            -snp-stats: Calculate and output per-SNP statistics.  This implies that no SNP filtering opt-
                                        ions are used.
                -snp-stats-columns <a>: Comma-seperated list of extra columns to output in the snp-wise statistics file.
                                        The standard columns are: SNPID, RSID, position, minor_allele, major_allele, MAF,
                                        HWE, missing, information. Your choices here are old_information, jonathans_info-
                                        rmation, mach_r2, and entropy.  Defaults to "".
                   -snp-stats-file <a>: Override the auto-generated path(s) of the snp-stats file to use when outputting
                                        snp-wise statistics.  (By default, the paths are formed by adding ".snp-stats" t-
                                        o the input gen filename(s).)  The '#' character can also be used here to specif-
                                        y one output file per chromosome.

SNP filtering options:
                              -hwe <a>: Filter out SNPs with -log10( HWE p-value ) greater than or equal to the value sp-
                                        ecified.
                         -info <a> <b>: Filter out SNPs with Fisher information lying outside the given range.
                          -maf <a> <b>: Filter out SNPs whose minor allele frequency lies outside the interval [a,b].
                -missing-call-rate <a>: Filter out SNPs with missing call rate greater than or equal to the value specif-
                                        ied.
                 -snp-interval <a> <b>: Filter out SNPs with position outside the interval [a,b].
                 -snp-missing-rate <a>: Filter out SNPs with missing data rate greater than or equal to the value specif-
                                        ied.

Sample filtering options:
               -heterozygosity <a> <b>: Filter out samples with heterozygosity outside the inteval [a,b].
              -sample-missing-rate <a>: Filter out samples with missing data rate greater than the value specified.

Inclusion / exclusion list options:
               -write-sample-excl-list: Do not apply sample filters directly.  Instead, write a file containing a list o-
                                        f the ids  of individuals which would be filtered out by the filter.
      -write-sample-excl-list-file <a>: Override default name of the file to use in -write-sample-excl-list
              -write-snp-excl-list <a>: Don't output a new genotypes file; instead, write files containing the positions
                                        of SNPs that are filtered out.
         -write-snp-excl-list-file <a>: Override the auto-generated path(s) of the file to use when outputting the posit-
                                        ions of filtered out SNPs.  (By default, the paths are formed by adding ".snp-ex-
                                        cl-list" to the input gen filename(s).)  If used, this option must appear as man-
                                        y times as the -g option.  If the corresponding occurence of -g uses a '#' wildc-
                                        ard character, the '#' character can also be used here to specify numbered outpu-
                                        t files corresponding to the input files.

Relatedness options:
                      -relatedness <a>: Compute relatedness matrices pairwise for all samples.  (This can take a long ti-
                                        me).
  -relatedness-alternative <a> <b> <c>: Set the probability of IBD0, IBD1, and IBD2 in the alternative model for the Bay-
                                        es factor. These must lie between 0 and 1 and sum to 1.  Defaults to "0, 0, 1".
              -relatedness-epsilon <a>: Set the probability of genotyping error at a SNP. This is used to make the model
                                        tolerant to genotyping errors.  Defaults to "0.001".
         -relatedness-null <a> <b> <c>: Set the probability of IBD0, IBD1, and IBD2 in the null model for the Bayes fact-
                                        or. These must lie between 0 and 1 and sum to 1.  Defaults to "1, 0, 0".
       -relatedness-sample-columns <a>: Choose ranges of samples to compute relatedness for. This option should be a com-
                                        ma-separated list of ranges of the form a-b, meaning use all rows between a and 
                                        b inclusive. A single value means a 1-element range; a range of the form a- or -a
                                        means all samples up to and including the given one.  Defaults to "0-".
          -relatedness-sample-rows <a>: Choose ranges of samples to compute relatedness for. This option should be a com-
                                        ma-separated list of ranges of the form a-b, meaning use all rows between a and 
                                        b inclusive. A single value means a 1-element range; a range of the form a- or -a
                                        means all samples up to and including the given one.  Defaults to "0-".

Association test options:
                       -covariates <a>: Specify a comma-separated list of covariates to use in the association test.  De-
                                        faults to "".  Currently this is ignored.
                               -ot <a>: Override the default name of the association test output file.  Defaults to "qct-
                                        ool.test".
                             -test <a>: Perform an association test on the given phenotype.  This must be a binary phenotype.

Other options:
                                -T <a>: Specify the number of worker threads to use in computationally intensive tasks. 
                                        Defaults to "0".
                                -force: Ignore warnings and proceed with requested action.
                              -log <a>: Override the default path of the log file written by qctool.  By default, this i-
                                        s qctool.log.  Defaults to "qctool.log".
                             -plot <a>: Path of file to produce plots in.


USAGE: gen-convert
==================

gen-convert converts one or more input gen files into one or more output gen files (possibly in a different gen file format).
Typical usage is:

> gen-convert -g (input files) -og (output files)

The -g option must appear the same number of times as the -og option.
As with qctool, a '#' wildcard character may appear in the input file and matches a number from 1 to 100.
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
