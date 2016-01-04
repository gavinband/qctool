
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <sstream>
#include <vector>
#include <string>
#include <set>
#include "genfile/GenomePosition.hpp"
#include "genfile/VCFFormatSNPDataSource.hpp"
#include "genfile/vcf/CallReader.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/get_set.hpp"
#include "test_case.hpp"
#include "CVCFT/CVCFT.cpp"

namespace data {
	// This example is from http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
	std::string const metadata4_0 =
		"##fileformat=VCFv4.0\n"
		"##FORMAT=<ID=GT,Number=.,Type=String,Description=\"Genotype\">\n"
		"##FORMAT=<ID=GL1,Number=3,Type=Float,Description=\"Genotype Likelihoods\">\n"
		"##FORMAT=<ID=GL3,Number=3,Type=Float,Description=\"Genotype Likelihoods\">\n"
		"##FORMAT=<ID=GL4,Number=4,Type=Float,Description=\"Genotype Likelihoods\">\n"
	;

	std::string const metadata4_1 =
		"##fileformat=VCFv4.1\n"
		"##FORMAT=<ID=GT,Number=.,Type=String,Description=\"Genotype\">\n"
		"##FORMAT=<ID=GL1,Number=3,Type=Float,Description=\"Genotype Likelihoods\">\n"
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

namespace impl {
	struct VectorSetter: public genfile::VariantDataReader::PerSampleSetter
	{
		VectorSetter( std::vector< double >& values ):
			m_values( values ),
			m_sample_i( 0 )
		{}
			
		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			// do nothing.
			m_sample_i = 0 ;
			m_values.clear() ;
		}
		bool set_sample( std::size_t i ) {
			m_sample_i = i ;
			return true ;
		}
		void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
			// do nothing
		}
		void set_value( std::size_t, MissingValue const value ) {
			m_values.push_back( 0.0 ) ;
		}
		void set_value( std::size_t, std::string& value ) {
			assert(0) ;
		}
		void set_value( std::size_t, Integer const value ) {
			assert(0) ;
		}
		void set_value( std::size_t, double const value ) {
			m_values.push_back( value ) ;
		}
		void finalise() {}
			
	public:
		std::vector< double >& m_values ;
		std::size_t m_sample_i ;
	} ;
}

void perform_test_against_cvcft(
	std::string const& data,
	std::size_t const number_of_lines,
	std::string const& key,
	std::size_t const number_of_probs_per_snp
) {
	using namespace data ;
	{
		std::vector< genfile::VariantIdentifyingData > cvcft_snps ;
		std::vector< std::vector< double > > cvcft_probs ;
		std::cerr << "CVCFT:\n" ;
		{
			std::auto_ptr< std::istream > istr( new std::istringstream( data ) ) ;
			CVCFT data( *istr, key.c_str() ) ;

			while( data.ReadNext() ) {
				genfile::GenomePosition position ;
				{
					std::string chr ;
					int pos ;
					data.GetChrPos( chr, pos ) ;
					position = genfile::GenomePosition( chr, pos ) ;
				}
				genfile::VariantIdentifyingData snp(
					std::string( data.GetRSID() ),
					position,
					std::string( 1, data.GetRef() ),
					std::string( 1, data.GetAlt() )
				) ;
				std::vector< double > prob ;
				data.GetProb( prob ) ;
				std::cerr << snp << ": " << key << ": " << prob.size() << " probs.\n" ;
				
				cvcft_snps.push_back( snp ) ;
				cvcft_probs.push_back( prob ) ;
			}

			TEST_ASSERT( cvcft_snps.size() == number_of_lines ) ;
		}

		std::vector< genfile::VariantIdentifyingData > my_snps ;
		std::vector< std::vector< double > > my_probs ;

		using namespace genfile::vcf ;

		std::cerr << "genfile::VCFFormatSNPDataSource:\n" ;
		{
			std::auto_ptr< std::istream > istr( new std::istringstream( data ) ) ;

			genfile::VCFFormatSNPDataSource source( istr ) ;
			TEST_ASSERT( !source.total_number_of_snps() || ( source.total_number_of_snps().get() == number_of_lines )) ;
			genfile::VariantIdentifyingData snp ;
			while( source.get_snp_identifying_data( &snp )) {
				std::vector< double > probs ;
				impl::VectorSetter setter( probs ) ;
				source.read_variant_data()->get( key, setter ) ;
				std::cerr << snp << ": " << key << ": " << probs.size() << " probs.\n" ;
				
				TEST_ASSERT( probs.size() == 0 || probs.size() == number_of_probs_per_snp ) ;
				my_snps.push_back( snp ) ;
				my_probs.push_back( probs ) ;
			}
		}
		
		TEST_ASSERT( my_snps.size() == my_probs.size() ) ;
		TEST_ASSERT( my_probs.size() == cvcft_probs.size() ) ;
		
		for( std::size_t i = 0; i < my_snps.size(); ++i ) {
			TEST_ASSERT( my_probs[i] == cvcft_probs[i] ) ;
			BOOST_CHECK_EQUAL( my_snps[i].get_rsid(), cvcft_snps[i].get_rsid() ) ;
			BOOST_CHECK_EQUAL( my_snps[i].get_position(), cvcft_snps[i].get_position() ) ;
			BOOST_CHECK_EQUAL( my_snps[i].get_allele(0), cvcft_snps[i].get_allele(0) ) ;
			BOOST_CHECK_EQUAL( my_snps[i].get_allele(1), cvcft_snps[i].get_allele(1) ) ;
		}
	}
}

AUTO_TEST_CASE( test_against_cvcft1 ) {
	perform_test_against_cvcft( data::metadata4_0 + data::data1, 7, "GL1", 4*3 ) ;
}
AUTO_TEST_CASE( test_against_cvcft2 ) {
	perform_test_against_cvcft( data::metadata4_0 + data::data1, 7, "GL3", 4*3 ) ;
}
AUTO_TEST_CASE( test_against_cvcft3 ) {
	perform_test_against_cvcft( data::metadata4_0 + data::data2, 5, "GL1", 4*3 ) ;
}
AUTO_TEST_CASE( test_against_cvcft4 ) {
	perform_test_against_cvcft( data::metadata4_0 + data::data2, 5, "GL3", 4*3 ) ;
}
AUTO_TEST_CASE( test_against_cvcft5 ) {
	perform_test_against_cvcft( data::metadata4_1 + data::data1, 7, "GL1", 4*3 ) ;
}
AUTO_TEST_CASE( test_against_cvcft6 ) {
	perform_test_against_cvcft( data::metadata4_1 + data::data1, 7, "GL3", 4*3 ) ;
}
AUTO_TEST_CASE( test_against_cvcft7 ) {
	perform_test_against_cvcft( data::metadata4_1 + data::data2, 5, "GL1", 4*3 ) ;
}
AUTO_TEST_CASE( test_against_cvcft8 ) {
	perform_test_against_cvcft( data::metadata4_1 + data::data2, 5, "GL3", 4*3 ) ;
}

