genfile, v1.0
-------------

genfile is a library for dealing with files in the GEN, BGEN, and BZGEN (BGEN with compressed SNP block) formats.

Also compared is bgen.gz, the zipped-up gen file format.


FORMAT COMPARISON
-----------------

All comparisons are made using gen-select revision 55 on my macbook pro (2.66GHz intel core duo with 7200rpm seagate drive).
The data used consists of a file with 34499 snps and probability data for 1504 individuals.

0. File sizes
source			size in bytes		size in millions of bytes (Mb)
------			-------------		------------------------------
gen				416,973,919			417Mb
bgen			312,431,645			312Mb
bzgen			52,297,230			52Mb
bgen.gz			43,532,732			44Mb

1. Convert a plain GEN file into a target file:
source		target		time in s
------		------		---------
gen			gen			153.9
gen	 		bgen		81.3
gen			bzgen		103.1
gen			bgen.gz		102.8
bgen		gen			91.9
bgen		bgen		19.9
bzgen		bgen		21.4
bgen		bzgen		41.0
bgen		bgen.gz		41.5

2. Time to read in a file using gen-select.  This does nothing except gen-select's standard processing.
source					time in s
------					---------
gen						70.4
bgen					10.2
bzgen					11.5
bgen.gz					12.5

3. Read times using benchmark-io (timing reads of 10000 SNPs from the given file)
source					time in s
------					---------
gen						20.0
bgen					2.3
bzgen					2.6
bgen.gz					2.8


CONCLUSIONS
-----------
1. gen format is slow both for reading and writing
2. bgen format is the fastest at both
3. For reading, bgen, bzgen and bgen.gz are comparable in speed with bzgen in the middle between bgen and bgen.gz
4. For writing, bzgen and bgen.gz are comparable, and may be roughly 3 times slower than bgen format.
5. For space, bzgen and bgen.gz use about 1/6th of the space of bgen and 1/8th the space of GEN.

RECOMMENDATION
--------------
* Use bzgen format unless space is not an issue, in which case use bgen.
