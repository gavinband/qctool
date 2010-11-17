#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "genfile/CategoricalCohortIndividualSource.hpp"
#include "genfile/SampleFilteringCohortIndividualSource.hpp"
#include <cassert>

#define TEST_ASSERT( thing ) assert( thing )

struct SampleFilteringCohortIndividualSourceTester
{
	SampleFilteringCohortIndividualSourceTester() {
		basic_test() ;
	}
	
	void basic_test() {
		std::cerr << "basic_test()\n" ;
		std::string sample_data = "ID_1 ID_2 missing C1 C2 P1 P2\n"
		"0 0 0 D C P B\n"
		"S1 S1 0.0 alpha 10.9 50 case\n"
		"S2 S2 0.0 alpha 10.9 50 case\n"
		"S3 S3 0.0 alpha 10.9 50 case\n"
		"S4 S4 0.0 alpha 10.9 50 case\n"
		"S5 S5 0.0 alpha 10.9 50 case\n"
		"S6 S6 0.0 alpha 10.9 50 case\n"
		"S7 S7 0.0 alpha 10.9 50 case\n"
		"S8 S8 0.0 alpha 10.9 50 case\n"
		"S9 S9 0.0 alpha 10.9 50 case\n"
		"S10 S10 0.0 alpha 10.9 50 case\n" ;

		std::vector< std::string > sample_ids ;
		sample_ids.push_back( "S1" ) ;
		sample_ids.push_back( "S2" ) ;
		sample_ids.push_back( "S3" ) ;
		sample_ids.push_back( "S4" ) ;
		sample_ids.push_back( "S5" ) ;
		sample_ids.push_back( "S6" ) ;
		sample_ids.push_back( "S7" ) ;
		sample_ids.push_back( "S8" ) ;
		sample_ids.push_back( "S9" ) ;
		sample_ids.push_back( "S10" ) ;

		std::vector< std::set< std::size_t > > test_exclusions ;
		test_exclusions.push_back( std::set< std::size_t >() ) ;
		test_exclusions.push_back( std::set< std::size_t >() ) ;
		test_exclusions.back().insert( 5 ) ;
		test_exclusions.push_back( std::set< std::size_t >() ) ;
		test_exclusions.back().insert( 5 ) ;
		test_exclusions.back().insert( 2 ) ;
		test_exclusions.push_back( std::set< std::size_t >() ) ;
		test_exclusions.back().insert( 5 ) ;
		test_exclusions.back().insert( 2 ) ;
		test_exclusions.back().insert( 9 ) ;
		test_exclusions.back().insert( 3 ) ;
		test_exclusions.push_back( std::set< std::size_t >() ) ;
		test_exclusions.back().insert( 0 ) ;
		test_exclusions.back().insert( 2 ) ;
		test_exclusions.back().insert( 4 ) ;
		test_exclusions.back().insert( 6 ) ;
		test_exclusions.back().insert( 8 ) ;

		for( std::size_t i = 0; i < test_exclusions.size(); ++i ) {
			genfile::CohortIndividualSource::UniquePtr source = genfile::CohortIndividualSource::create(
				"data://" + sample_data,
				"categorical"
			) ;
			genfile::CohortIndividualSource::ColumnSpec column_spec = source->get_column_spec() ;
			std::size_t N = source->get_number_of_individuals() ;
			
			genfile::SampleFilteringCohortIndividualSource::UniquePtr filtering_source = genfile::SampleFilteringCohortIndividualSource::create(
				source,
				test_exclusions[i]
			) ;

			TEST_ASSERT( filtering_source->get_column_spec() == column_spec ) ;
			TEST_ASSERT( filtering_source->get_number_of_individuals() == N - test_exclusions[i].size() ) ;
			
			for( std::size_t j = 0, k = 0; j < filtering_source->get_number_of_individuals(); ++j, ++k ) {
				while( test_exclusions[i].find( k ) != test_exclusions[i].end() ) {
					++k ;
				}
				TEST_ASSERT( filtering_source->get_entry( j, "id_1" ).as< std::string >() == sample_ids[k] ) ;
			}
			std::cerr << "basic_test(): success.\n" ;
		}
	}
	
} ;

int main() {
	SampleFilteringCohortIndividualSourceTester tester ;
}
