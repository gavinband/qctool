# qctool - README
Copyright 2009-2021 Gavin Band

This repository contains the source code for QCTOOL and a number of other command-line programs
that manipulate gwas datasets and other genomic data. The repository is hosted [at code.enkre.net](https://code.enkre.net/qctool) and is also mirrored [on github](https://github.com/gavinband/qctool).  The programs included are:

* [QCTOOL](http://www.well.ox.ac.uk/~gav/qctool), a program for manipulation and quality control of gwas datasets and other genome-wide data.  Read [the QCTOOL documentation](http://www.well.ox.ac.uk/~gav/qctool).
* [Inthinnerator](http://www.well.ox.ac.uk/~gav/inthinnerator), a program for thinning genetic variants using a recombination map.  Read [the inthinnerator documentation](http://www.well.ox.ac.uk/~gav/inthinnerator).
* [HPTEST](http://www.well.ox.ac.uk/~gav/hptest), a program for testing for association between host and parasite genotypes. Read the [HPTEST documentation](http://www.well.ox.ac.uk/~gav/hptest) for more details.
* [LDBIRD](http://www.well.ox.ac.uk/~gav/ldbird), a program for computing LD metrics between all pairs of variants in a genomic region or across whole genomes.  Read the [LDBIRD documentation](http://www.well.ox.ac.uk/~gav/ldbird) for more details.

## Obtaining QCTOOL ##

Basic instructions for obtaining and/or compiling this package can be found on the
[QCTOOL download page](http://www.well.ox.ac.uk/~gav/qctool/documentation/download.html).
For further information, see the wiki page on [compiling QCTOOL](/wiki/Compiling%20QCTOOL).

## License ##

The QCTOOL package is copyright Gavin Band released under the [Boost software license](LICENSE_1_0.txt). QCTOOL uses
code from 3rd party packages described below, that use their own licenses. These include the BSD License and the GNU
General Public License. See the [respective folders for these packages](/dir?ci=tip&name=3rd_party) for details.

## Acknowledgements ##

QCTOOL was written by Gavin Band, based on an initial design by Gavin Band and Jonathan Marchini.

Inthinnerator, HPTEST, LDBIRD and all other included software was written by Gavin Band.

QCTOOL gratefully makes use of functionality from the following projects:

- zlib (<http://zlib.net/>)
- the boost C++ libraries (<http://www.boost.org>)
- sqlite3 (<http://www.sqlite.org>),
- the Eigen library (<http://eigen.tuxfamily.org>)
- the SNPHWE code from Wigginton et al (<http://www.sph.umich.edu/csg/abecasis/Exact>)
- the Zstandard library(<http://www.zstd.net>).

Quang Si Le contributed an implementation of the VCF format used for testing.
