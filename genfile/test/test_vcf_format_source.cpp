#include <sstream>
#include <vector>
#include <string>
#include <set>
#include "genfile/GenomePosition.hpp"
#include "genfile/VCFFormatSNPDataSource.hpp"
#include "genfile/vcf/CallReader.hpp"
#include "genfile/Error.hpp"
#include "test_case.hpp"
#include "CVCFT/CVCFT.cpp"

namespace data {
	// This example is from http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
	std::string const metadata4_0 =
		"##fileformat=VCFv4.0\n"
		"##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n"
		"##FORMAT=<ID=GL1,Number=3,Type=Float,Description=\"Genotype Likelihoods\">\n"
		"##FORMAT=<ID=GL3,Number=3,Type=Float,Description=\"Genotype Likelihoods\">\n"
		"##FORMAT=<ID=GL4,Number=4,Type=Float,Description=\"Genotype Likelihoods\">\n"
	;

	std::string const metadata4_1 =
		"##fileformat=VCFv4.1\n"
		"##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n"
		"##FORMAT=<ID=GL1,Number=.,Type=Float,Description=\"Genotype Likelihoods\">\n"
		"##FORMAT=<ID=GL3,Number=G,Type=Float,Description=\"Genotype Likelihoods\">\n"
		"##FORMAT=<ID=GL4,Number=4,Type=Float,Description=\"Genotype Likelihoods\">\n"
	;

	std::string const data1 = 
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND1\tIND2\tIND3\tIND4\n"
		"1\t1\trs1\tA\tG\t1\tPASS\t\tGT:GL1\t1|1:0.1,0.2,0.3\t0|0:0.5,0.6,0.7\t0|1:1.6627,1.7627,1.8627\t0/0:1.22e04,1.32e04,1.42e04\n"
		"1\t2\trs2\tC\tT\t1\tPASS\t\tGT:GL1\t1|1:1.1,1.2,1.3\t0|0:1.5,1.6,1.7\t0|1:2.6627,2.7627,2.8627\t0/0:2.22e05,2.32e04,2.42e04\n"
		"1\t3\trs3\tA\tC\t1\tPASS\t\tGT:GL1\t1|1:2.1,2.2,2.3\t0|0:2.5,2.6,2.7\t0|1:3.6627,3.7627,3.8627\t0/0:3.22e05,3.32e04,3.42e04\n"
		"2\t4\trs4\tT\tA\t1\tPASS\t\tGT:GL1:GL3\t1|1:2.1,2.2,2.3:1.0,1.1,1.2\t0|0:2.5,2.6,2.7:2.0,2.1,2.2\t0|1:3.6627,3.7627,3.8627:3.0,3.1,3.2\t0/0:3.22e05,3.32e05,3.42e05:4.0,4.1,4.2\n"
		"2\t4\trs4\tT\tA\t1\tPASS\t\tGT:GL3:GL1\t1|1:1.0,1.1,1.2:2.1,2.2,2.3\t0|0:2.0,2.1,2.2:2.5,2.6,2.7\t0|1:3.0,3.1,3.2:3.6627,3.7627,3.8627\t0/0:4.0,4.1,4.2:3.22e05,3.32e05,3.42e05\n"
		"2\t4\trs4\tT\tA\t1\tPASS\t\tGT:GL1:GL4\t1|1:2.1,2.2,2.3:1.0,1.1,1.2,0.0001\t0|0:2.5,2.6,2.7:2.0,2.1,2.2,0.0002\t0|1:3.6627,3.7627,3.8627:3.0,3.1,3.2,0.0003\t0/0:3.22e05,3.32e05,3.42e05:4.0,4.1,4.2,0.0004\n"
		"2\t4\trs4\tT\tA\t1\tPASS\t\tGT:GL4:GL1\t1|1:1.0,1.1,1.2,0.1001:2.1,2.2,2.3\t0|0:2.0,2.1,2.2,0.1002:2.5,2.6,2.7\t0|1:3.0,3.1,3.2,0.1003:3.6627,3.7627,3.8627\t0/0:4.0,4.1,4.2,0.1004:3.22e05,3.32e05,3.42e05\n"
	;

	std::string const data2 = 
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND1\tIND2\tIND3\tIND4\n"
		"1\t1\trs1\tA\tG\t1\tPASS\t\tGT:GL1\t1|1\t0|0\t0|1\t0/0\n"
		"1\t2\trs2\tC\tT\t1\tPASS\t\tGT:GL1\t1|1\t0|0\t0|1\t0/0\n"
		"1\t3\trs3\tA\tC\t1\tPASS\t\tGT:GL1\t1|1\t0|0\t0|1\t0/0\n"
		"2\t4\trs4\tT\tA\t1\tPASS\t\tGT:GL1\t1|1\t0|0\t0|1\t0/0:3.22e05,3.32e05,3.42e05\n"
		"2\t4\trs4\tT\tA\t1\tPASS\t\tGT:GL1\t1|1\t0|0\t0|1\t0/0:4.0,4.1,4.2\n"
	;
}

void perform_test_against_cvcft(
	std::string const& data,
	std::size_t const number_of_lines,
	std::string const& key,
	std::size_t const number_of_probs_per_snp
) {
	using namespace data ;
	{
		std::vector< genfile::SNPIdentifyingData > cvcft_snps ;
		std::vector< std::vector< double > > cvcft_probs ;
		{
			std::auto_ptr< std::istream > istr( new std::istringstream( data ) ) ;
			CVCFT data( *istr, key.c_str() ) ;

			while( data.ReadNext() ) {
				genfile::SNPIdentifyingData snp ;
				snp.rsid() = data.GetRSID() ;
				snp.first_allele() = data.GetRef() ;
				snp.second_allele() = data.GetAlt() ;
				{
					std::string chr ;
					int pos ;
					data.GetChrPos( chr, pos ) ;
					snp.position() = genfile::GenomePosition( chr, pos ) ;
				}
				std::vector< double > prob ;
				data.GetProb( prob ) ;
				
				cvcft_snps.push_back( snp ) ;
				cvcft_probs.push_back( prob ) ;
			}

			TEST_ASSERT( cvcft_snps.size() == number_of_lines ) ;
		}

		std::vector< genfile::SNPIdentifyingData > my_snps ;
		std::vector< std::vector< double > > my_probs ;

		using namespace genfile::vcf ;

		{
			std::auto_ptr< std::istream > istr( new std::istringstream( data ) ) ;

			genfile::VCFFormatSNPDataSource source( istr, key ) ;
			TEST_ASSERT( source.total_number_of_snps() == number_of_lines ) ;
			genfile::SNPIdentifyingData snp ;
			while( source.get_snp_identifying_data( snp )) {
				std::vector< double > probs ;
				source.read_snp_probability_data( genfile::set_genotypes( probs )) ;
				TEST_ASSERT( probs.size() == 0 || probs.size() == number_of_probs_per_snp ) ;
				my_snps.push_back( snp ) ;
				my_probs.push_back( probs ) ;
			}
			TEST_ASSERT( my_snps.size() == source.total_number_of_snps() ) ;
			TEST_ASSERT( my_probs.size() == source.total_number_of_snps() ) ;
		}
		
		TEST_ASSERT( my_snps.size() == my_probs.size() ) ;
		TEST_ASSERT( my_probs.size() == cvcft_probs.size() ) ;
		
		for( std::size_t i = 0; i < my_snps.size(); ++i ) {
			TEST_ASSERT( my_probs[i] == cvcft_probs[i] ) ;
			TEST_ASSERT( my_snps[i].get_rsid() == cvcft_snps[i].get_rsid() ) ;
			TEST_ASSERT( my_snps[i].get_position() == cvcft_snps[i].get_position() ) ;
			TEST_ASSERT( my_snps[i].get_first_allele() == cvcft_snps[i].get_first_allele() ) ;
			TEST_ASSERT( my_snps[i].get_second_allele() == cvcft_snps[i].get_second_allele() ) ;
		}
	}
}

AUTO_TEST_CASE( test_against_cvcft ) {
	std::cerr << "test_against_cvcft()..." ;
	perform_test_against_cvcft( data::metadata4_0 + data::data1, 7, "GL1", 4*3 ) ;
	perform_test_against_cvcft( data::metadata4_0 + data::data1, 7, "GL3", 4*3 ) ;
	perform_test_against_cvcft( data::metadata4_0 + data::data2, 5, "GL1", 4*3 ) ;
	perform_test_against_cvcft( data::metadata4_0 + data::data2, 5, "GL3", 4*3 ) ;

	perform_test_against_cvcft( data::metadata4_1 + data::data1, 7, "GL1", 4*3 ) ;
	perform_test_against_cvcft( data::metadata4_1 + data::data1, 7, "GL3", 4*3 ) ;
	perform_test_against_cvcft( data::metadata4_1 + data::data2, 5, "GL1", 4*3 ) ;
	perform_test_against_cvcft( data::metadata4_1 + data::data2, 5, "GL3", 4*3 ) ;
	std::cerr << "ok.\n" ;
}

AUTO_TEST_MAIN {
	test_against_cvcft() ;
}
