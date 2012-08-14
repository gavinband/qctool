
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <sstream>
#include "test_case.hpp"
#include "genfile/SNPIdentifyingData2.hpp"
#include "genfile/GenomePosition.hpp"

AUTO_TEST_CASE( test_output_streaming ) {
	for( char c1 = 33; c1 <= 126; ++c1 ) {
		std::string rsid = std::string( 1, c1 ) ;
		genfile::SNPIdentifyingData2 snp( rsid, genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
		std::ostringstream ostr ;
		ostr << snp ;
		TEST_ASSERT( ostr.str() == rsid + " NA 0 A G" ) ;
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
				TEST_ASSERT( ostr.str() == rsid + " NA 0 A G" ) ;
			}
			else {
				TEST_ASSERT( ostr.str() == rsid + " [" + alternate_id + "] NA 0 A G" ) ;
			}
		}
	}
}

AUTO_TEST_CASE( test_alternate_ids ) {
	for( std::size_t n = 0; n < 10; ++n ) {
		for( char c1 = 33; c1 <= 126; ++c1 ) {
			std::string rsid = std::string( 1, c1 ) ;
			genfile::SNPIdentifyingData2 snp( rsid, genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
			for( std::size_t i = 0; i < n; ++i ) {
				snp.add_identifier( std::string( i + 2, c1 ) ) ;
			}
			std::vector< genfile::string_utils::slice > ids = snp.get_identifiers() ;
			TEST_ASSERT( ids.size() == n ) ;
			for( std::size_t i = 0; i < n; ++i ) {
				TEST_ASSERT( ids[i] == std::string( i + 2, c1 ) ) ;
			}
		}
	}
}

