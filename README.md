# qctool - README
Copyright 2009-2016 Gavin Band, University of Oxford

This package comprises the QCTOOL program for quality control of a set of SNP marker data,
together with some other utility programs.  For full details see <http://www.well.ox.ac.uk/~gav/qctool>.

## INSTALLATION

### PREREQUISITES

QCTOOL requires the following external libraries which must be installed (with headers before compilation):

- zlib (http://zlib.net/)

QCTOOL also depends on the boost libraries (<http://www.boost.org>), sqlite3 (<http://www.sqlite.org>),
the Eigen library (<http://eigen.tuxfamily.org>) as well as the SNPHWE code from Wigginton et al
(http://www.sph.umich.edu/csg/abecasis/Exact).  All of these are included in the qctool repository and don't need
to installed beforehand.

### Steps to compile
Compilation involves two steps: configuration and build.
To configure the package, from the top-level directory run
`./waf-1.5.18 configure`
Or:
`./waf-1.5.18 configure --prefix path/to/install/location/`
to specify an installation location other than /usr/local.

To build the package, run
`./waf-1.5.18`
Then main qctool executable can be found under build/release/.  For debugging purposes, a debug build is also built under build/default/.

To install the executable, run
`./waf-1.5.18 install`