
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <sstream>
#include <vector>
#include <boost/bind.hpp>
#include "test_case.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPIdentifyingData2.hpp"
#include "genfile/GenomePosition.hpp"

AUTO_TEST_CASE( test_output_streaming ) {
	for( char c1 = 33; c1 <= 126; ++c1 ) {
		std::string rsid = std::string( 1, c1 ) ;
		genfile::SNPIdentifyingData2 snp( rsid, genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
		std::ostringstream ostr ;
		ostr << snp ;
		BOOST_CHECK_EQUAL( ostr.str(), rsid + " NA 0 A G" ) ;
	}

	for( char c1 = 33; c1 <= 126; ++c1 ) {
		for( char c2 = 33; c2 <= 126; ++c2 ) {
			std::string rsid = std::string( 1, c1 ) ;
			std::string alternate_id = std::string( 1, c2 ) ;
			genfile::SNPIdentifyingData2 snp( rsid, genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
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
			genfile::SNPIdentifyingData2 snp( rsid, genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
			snp.add_identifier( rsid ) ; // should be a no-op.
			for( std::size_t i = 0; i < n; ++i ) {
				snp.add_identifier( std::string( i + 2, c1 ) ) ;
			}
			std::vector< genfile::string_utils::slice > ids = snp.get_alternative_identifiers() ;
			BOOST_CHECK_EQUAL( ids.size(), n ) ;
			for( std::size_t i = 0; i < n; ++i ) {
				BOOST_CHECK_EQUAL( ids[i], std::string( i + 2, c1 ) ) ;
			}
		}
	}
}

AUTO_TEST_CASE( test_snp_values ) {
	genfile::SNPIdentifyingData2 snp( "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 1000 ), "A", "G" ) ;
	BOOST_CHECK_EQUAL( snp.get_rsid(), "RSID_1" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), genfile::GenomePosition( genfile::Chromosome(), 1000 ) ) ;
	BOOST_CHECK_EQUAL( snp.get_first_allele(), "A" ) ;
	BOOST_CHECK_EQUAL( snp.get_second_allele(), "G" ) ;
	BOOST_CHECK( snp.get_alternative_identifiers() == std::vector< genfile::string_utils::slice >() ) ;
}

AUTO_TEST_CASE( test_snp_copying ) {
	genfile::SNPIdentifyingData2 snp1( "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 1000 ), "A", "G" ) ;
	snp1.add_identifier( std::string( "SNP1" ) ) ;

	genfile::SNPIdentifyingData2 snp2 = snp1 ;
	BOOST_CHECK_EQUAL( snp1, snp2 ) ;

	genfile::SNPIdentifyingData snp3( "SNP1", "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 1000 ), "A", "G" ) ;
	genfile::SNPIdentifyingData2 snp4 = snp3 ;
	BOOST_CHECK_EQUAL( snp4, snp2 ) ;
}

AUTO_TEST_CASE( test_alternate_ids ) {
	std::string rsid ;
	{
		genfile::SNPIdentifyingData2 snp( "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 0 ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 0 ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 1 ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "SNPID_1" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 1 ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "SNPID_1" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 1 ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "SNPID_1" ) ;
		snp.add_identifier( std::string( "AnotherID" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[1], "AnotherID" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[1], "AnotherID" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[1], "AnotherID" ) ;
	}
}

AUTO_TEST_CASE( test_alternate_ids_2 ) {
	std::string rsid ;
	std::vector< std::string > ids ;
	{
		genfile::SNPIdentifyingData2 snp( "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
		ids.clear() ;
		snp.get_alternative_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 0 ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		ids.clear() ;
		snp.get_alternative_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 0 ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		ids.clear() ;
		snp.get_alternative_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 1 ) ;
		BOOST_CHECK_EQUAL( ids[0], "SNPID_1" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		ids.clear() ;
		snp.get_alternative_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 1 ) ;
		BOOST_CHECK_EQUAL( ids[0], "SNPID_1" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		ids.clear() ;
		snp.get_alternative_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 1 ) ;
		BOOST_CHECK_EQUAL( ids[0], "SNPID_1" ) ;
		snp.add_identifier( std::string( "AnotherID" ) ) ;
		ids.clear() ;
		snp.get_alternative_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 2 ) ;
		BOOST_CHECK_EQUAL( ids[0], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "AnotherID" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		ids.clear() ;
		snp.get_alternative_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 2 ) ;
		BOOST_CHECK_EQUAL( ids[0], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "AnotherID" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		ids.clear() ;
		snp.get_alternative_identifiers( boost::bind( &std::vector< std::string >::push_back, &ids, _1 ) ) ;
		BOOST_CHECK_EQUAL( ids.size(), 2 ) ;
		BOOST_CHECK_EQUAL( ids[0], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( ids[1], "AnotherID" ) ;
	}
}

AUTO_TEST_CASE( test_snp_data_setters ) {
	genfile::GenomePosition pos( genfile::Chromosome( "01" ), 0 ) ;
	genfile::SNPIdentifyingData2 snp( "RSID_1", pos, "A", "G" ) ;
	snp.add_identifier( std::string( "ID" ) ) ;
	snp.set_rsid( std::string( "RSID_1_changed" ) ) ;
	BOOST_CHECK_EQUAL( snp.get_rsid(), "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_first_allele(), "A" ) ;
	BOOST_CHECK_EQUAL( snp.get_second_allele(), "G" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
	BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 1 ) ;
	BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "ID" ) ;

	snp.set_first_allele( std::string( "GA" ) ) ;
	BOOST_CHECK_EQUAL( snp.get_rsid(), "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_first_allele(), "GA" ) ;
	BOOST_CHECK_EQUAL( snp.get_second_allele(), "G" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
	BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 1 ) ;
	BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "ID" ) ;

	snp.set_second_allele( std::string( "AG" ) ) ;
	BOOST_CHECK_EQUAL( snp.get_rsid(), "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_first_allele(), "GA" ) ;
	BOOST_CHECK_EQUAL( snp.get_second_allele(), "AG" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
	BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 1 ) ;
	BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "ID" ) ;

	pos.position() == 1000 ;
	snp.set_position( pos ) ;
	BOOST_CHECK_EQUAL( snp.get_rsid(), "RSID_1_changed" ) ;
	BOOST_CHECK_EQUAL( snp.get_first_allele(), "GA" ) ;
	BOOST_CHECK_EQUAL( snp.get_second_allele(), "AG" ) ;
	BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
	BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 1 ) ;
	BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "ID" ) ;
}

