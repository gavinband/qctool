
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <sstream>
#include <vector>
#include <boost/bind.hpp>
#include "test_case.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/GenomePosition.hpp"

AUTO_TEST_CASE( test_output_streaming ) {
	for( char c1 = 33; c1 <= 126; ++c1 ) {
		std::string rsid = std::string( 1, c1 ) ;
		genfile::VariantIdentifyingData snp( rsid, genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
		std::ostringstream ostr ;
		ostr << snp ;
		BOOST_CHECK_EQUAL( ostr.str(), rsid + " NA 0 A G" ) ;
	}

	for( char c1 = 33; c1 <= 126; ++c1 ) {
		for( char c2 = 33; c2 <= 126; ++c2 ) {
			std::string rsid = std::string( 1, c1 ) ;
			std::string alternate_id = std::string( 1, c2 ) ;
			genfile::VariantIdentifyingData snp( rsid, genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
			snp.add_identifier( alternate_id ) ;
			std::ostringstream ostr ;
			ostr << snp ;
			if( c1 == c2 ) {
				BOOST_CHECK_EQUAL( ostr.str(), rsid + " NA 0 A G" ) ;
			}
			else {
				BOOST_CHECK_EQUAL( ostr.str(), rsid + " [" + alternate_id + "] NA 0 A G" ) ;
			}
		}
	}
}

AUTO_TEST_CASE( test_alternate_ids2 ) {
	for( std::size_t n = 0; n < 10; ++n ) {
		for( char c1 = 33; c1 <= 126; ++c1 ) {
			std::string rsid = std::string( 1, c1 ) ;
			genfile::VariantIdentifyingData snp( rsid, genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
			snp.add_identifier( rsid ) ; // should be a no-op.
			for( std::size_t i = 0; i < n; ++i ) {
				snp.add_identifier( std::string( i + 2, c1 ) ) ;
			}
			std::vector< genfile::string_utils::slice > ids = snp.get_identifiers( 1, snp.number_of_identifiers() ) ;
			BOOST_CHECK_EQUAL( ids.size(), n ) ;
			for( std::size_t i = 0; i < n; ++i ) {
				BOOST_CHECK_EQUAL( ids[i], std::string( i + 2, c1 ) ) ;
			}
		}
	}
}

AUTO_TEST_CASE( test_snp_values ) {
	genfile::VariantIdentifyingData snp( "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 1000 ), "A", "G" ) ;
	BOOST_CHECK_EQUAL( snp.get_primary_id(), "RSID_1" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), genfile::GenomePosition( genfile::Chromosome(), 1000 ) ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(0), "A" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(1), "G" ) ;
	BOOST_CHECK( snp.get_identifiers() == std::vector< genfile::string_utils::slice >( 1, "RSID_1" ) ) ;
}

AUTO_TEST_CASE( test_snp_copying ) {
	genfile::VariantIdentifyingData snp1( "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 1000 ), "A", "G" ) ;
	snp1.add_identifier( std::string( "SNP1" ) ) ;

	genfile::VariantIdentifyingData snp2 = snp1 ;
	BOOST_CHECK_EQUAL( snp1, snp2 ) ;

	genfile::VariantIdentifyingData snp3( "SNP1", "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 1000 ), "A", "G" ) ;
	genfile::VariantIdentifyingData snp4 = snp3 ;
	BOOST_CHECK_EQUAL( snp4, snp2 ) ;
}

AUTO_TEST_CASE( test_alternate_ids_2 ) {
	std::string rsid ;
	std::vector< std::string > ids ;
	{
		genfile::VariantIdentifyingData snp( "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
		ids.clear() ;
		snp.get_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 1 ) ;
		BOOST_CHECK_EQUAL( ids[0], "RSID_1" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		ids.clear() ;
		snp.get_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 1 ) ;
		BOOST_CHECK_EQUAL( ids[0], "RSID_1" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		ids.clear() ;
		snp.get_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 2 ) ;
		BOOST_CHECK_EQUAL( ids[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "SNPID_1" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		ids.clear() ;
		snp.get_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 2 ) ;
		BOOST_CHECK_EQUAL( ids[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "SNPID_1" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		ids.clear() ;
		snp.get_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 2 ) ;
		BOOST_CHECK_EQUAL( ids[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "SNPID_1" ) ;
		snp.add_identifier( std::string( "AnotherID" ) ) ;
		ids.clear() ;
		snp.get_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 3 ) ;
		BOOST_CHECK_EQUAL( ids[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( ids[2], "AnotherID" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		ids.clear() ;
		snp.get_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 3 ) ;
		BOOST_CHECK_EQUAL( ids[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( ids[2], "AnotherID" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		ids.clear() ;
		snp.get_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 3 ) ;
		BOOST_CHECK_EQUAL( ids[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( ids[2], "AnotherID" ) ;
	}
}

AUTO_TEST_CASE( test_snp_data_setters ) {
	genfile::GenomePosition pos( genfile::Chromosome( "01" ), 0 ) ;
	genfile::VariantIdentifyingData snp( "RSID_1", pos, "A", "G" ) ;
	snp.add_identifier( std::string( "ID" ) ) ;
	snp.set_primary_id( std::string( "RSID_1_changed" ) ) ;
	BOOST_CHECK_EQUAL( snp.get_primary_id(), "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(0), "A" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(1), "G" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 2 ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "ID" ) ;

	snp.set_allele( 0, std::string( "GA" ) ) ;
	BOOST_CHECK_EQUAL( snp.get_primary_id(), "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(0), "GA" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(1), "G" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 2 ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "ID" ) ;

	snp.set_allele( 1, std::string( "AG" ) ) ;
	BOOST_CHECK_EQUAL( snp.get_primary_id(), "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(0), "GA" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(1), "AG" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 2 ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "ID" ) ;

	pos.position() = 1000 ;
	snp.set_position( pos ) ;
	BOOST_CHECK_EQUAL( snp.get_primary_id(), "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(0), "GA" ) ;
	BOOST_CHECK_EQUAL( snp.get_allele(1), "AG" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 2 ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "ID" ) ;
}

