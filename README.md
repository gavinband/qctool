# qctool - README
Copyright 2009-2016 Gavin Band, University of Oxford

This package comprises the QCTOOL program for quality control of a set of SNP marker data,
together with some other utility programs.  For full details see <http://www.well.ox.ac.uk/~gav/qctool>.

## INSTALLATION

### PREREQUISITES

QCTOOL requires the following external libraries which must be installed (with headers before compilation):

- zlib (http://zlib.net/)

Zlib is already installed on most systems.

QCTOOL also depends on the boost libraries (<http://www.boost.org>), sqlite3 (<http://www.sqlite.org>),
the Eigen library (<http://eigen.tuxfamily.org>) as well as the SNPHWE code from Wigginton et al
(http://www.sph.umich.edu/csg/abecasis/Exact).  All of these are included in the qctool repository and don't need
to installed beforehand.

## Compiling QCTOOL ##

Basic instructions for downloading and compiling QCTOOL can be found on the [QCTOOL download page](http://www.well.ox.ac.uk/~gav/qctool_v2/#download).  For further information, see the wiki page on [compiling QCTOOL](https://bitbucket.org/gavinband/qctool/wiki/Compiling%20QCTOOL).