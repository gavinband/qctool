#include <sstream>
#include <vector>
#include <string>
#include <set>
#include "genfile/VCFFormatSNPDataSource.hpp"
#include "genfile/Error.hpp"
#include "test_case.hpp"

namespace data {
	// This example is from http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
	std::string const ex1 =
		"##fileformat=VCFv4.1\n"
		"##fileDate=20090805\n"
		"##source=myImputationProgramV3.1\n"
		"##reference=file:/seq/references/1000GenomesPilot-NCBI36.fasta\n"
		"##contig=<name=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=Homo sapiens,taxonomy=x>\n"
		"##phasing=partial\n" ;
	
	std::string const ex2 =
		"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
		"##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"
		"##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">\n"
		"##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">\n" ;
	
	std::string const ex3 =
		"##FILTER=<ID=q10,Description=\"Quality below 10\">\n"
		"##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n" ;
		
	std::string const ex4 =
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
		"##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n"
		"##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Genotype Likelihoods\">\n" ;
		
	std::string const ex5 =
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n" ;

	std::string const ex6 = 
		"20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.\n"
		"20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3\n"
		"20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4\n"
		"20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2\n"
		"20\t1234567\tmicrosat1\tGTC\tG,GTCTC\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3\n" ;
}

AUTO_TEST_CASE( test_spec_example ) {
	using namespace data ;
	std::cerr << "test_spec_example()..." ;
	{
		std::auto_ptr< std::istream > istr( new std::istringstream( ex1 + ex2 + ex3 + ex4 + ex5 + ex6 ) ) ;
		genfile::VCFFormatSNPDataSource source( istr, "GL" ) ;
		TEST_ASSERT( source.total_number_of_snps() == 5 ) ;
	}
	{
		std::auto_ptr< std::istream > istr( new std::istringstream( ex1 + ex2 + ex3 + ex4 + ex5 ) ) ;
		genfile::VCFFormatSNPDataSource source( istr, "GL" ) ;
		TEST_ASSERT( source.total_number_of_snps() == 0 ) ;
	}
	std::cerr << "ok.\n" ;
}


AUTO_TEST_MAIN {
	test_spec_example() ;
}
