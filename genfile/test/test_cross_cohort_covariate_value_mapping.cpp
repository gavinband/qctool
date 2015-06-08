
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include "test_case.hpp"
#include "genfile/TraditionalStrictCohortIndividualSource.hpp"
#include "genfile/CrossCohortCovariateValueMapping.hpp"
#include "genfile/CategoricalCrossCohortCovariateValueMapping.hpp"
#include "genfile/ContinuousVariableCrossCohortCovariateValueMapping.hpp"
#include "genfile/LevelCountingCrossCohortCovariateValueMapping.hpp"
#include "genfile/NormalisingCrossCohortCovariateValueMapping.hpp"
#include "genfile/string_utils.hpp"
#include <cassert>

struct CategoricalValueMappingTester
{
	CategoricalValueMappingTester() {
		number_of_mapped_values_test() ;
		level_value_test() ;
		multiple_cohort_test() ;
	}
	
	// Test that with a single source we get the correct number of levels.
	void number_of_mapped_values_test() {
		std::cerr << "CategoricalValueMappingTester::number_of_mapped_values_test()\n" ;
		std::vector< std::string > test_data ;
		std::vector< std::size_t > expected_number_of_levels ;
		
		// 0 -entry tests
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D" ) ;
		expected_number_of_levels.push_back( 0 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 B" ) ;
		expected_number_of_levels.push_back( 0 ) ;

		// Discrete covariate tests
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		expected_number_of_levels.push_back( 1 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1\nS2 S2 0 2" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1\nS2 S2 0 3" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 3" ) ;
		expected_number_of_levels.push_back( 3 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 2" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 2\nS4 S4 0 1" ) ;
		expected_number_of_levels.push_back( 2 ) ;

		
		// binary phenotype tests
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 B\nS1 S1 0 1" ) ;
		expected_number_of_levels.push_back( 1 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 B\nS1 S1 0 1\nS2 S2 0 0" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 B\nS1 S1 0 1\nS2 S2 0 1" ) ;
		expected_number_of_levels.push_back( 1 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 B\nS1 S1 0 1\nS2 S2 0 1\nS3 S3 0 0" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 B\nS1 S1 0 1\nS2 S2 0 0\nS3 S3 0 0" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		
		// For each of the above, try again with a number of missing values added.
		{
			std::size_t const N = test_data.size() ;
			for( std::size_t i = 0; i < N; ++i ) {
				std::string data = test_data[i] ;
				for( std::size_t j = 1; j < 4; ++j ) {
					data += "\nSM" + genfile::string_utils::to_string( j ) + " SM" + genfile::string_utils::to_string( j ) + " 0 NA" ;
				}
				test_data.push_back( data ) ;	
				expected_number_of_levels.push_back( expected_number_of_levels[i] ) ;
			}
		}
		
		// Run the tests.
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "test #" << i << "\n" ;
			std::istringstream istr( test_data[i] ) ;
			genfile::TraditionalStrictCohortIndividualSource source( istr ) ;
			genfile::CategoricalCrossCohortCovariateValueMapping mapping( "column1" ) ;
			mapping.add_source( source ) ;
			TEST_ASSERT( mapping.get_number_of_distinct_mapped_values() == expected_number_of_levels[i] ) ;
		}
	}

	void level_value_test() {
		std::cerr << "CategoricalValueMappingTester::level_value_test()\n" ;
		
		std::vector< std::string > test_data ;
		std::vector< genfile::CohortIndividualSource::Entry > expecteds ;
		
		// In all cases the levels should be NA,1,2,3,4.
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 3\nS4 S4 0 4\nS5 S5 0 NA" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 NA\nS2 S2 0 4\nS3 S3 0 3\nS4 S4 0 2\nS5 S5 0 1" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 4\nS2 S2 0 1\nS3 S3 0 2\nS4 S4 0 3\nS5 S5 0 NA" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2" ) ;

		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 3\nS4 S4 0 4\nS5 S5 0 NA\nS6 S6 0 1" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 NA\nS2 S2 0 4\nS3 S3 0 3\nS4 S4 0 2\nS5 S5 0 1\nS6 S6 0 2" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 4\nS2 S2 0 1\nS3 S3 0 2\nS4 S4 0 3\nS5 S5 0 NA\nS6 S6 0 3" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2\nS6 S6 0 4" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2\nS6 S6 0 NA" ) ;

		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 3\nS4 S4 0 4\nS5 S5 0 NA\nS6 S6 0 1\nS7 S7 0 2" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 NA\nS2 S2 0 4\nS3 S3 0 3\nS4 S4 0 2\nS5 S5 0 1\nS6 S6 0 2\nS7 S7 0 1" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 4\nS2 S2 0 1\nS3 S3 0 2\nS4 S4 0 3\nS5 S5 0 NA\nS6 S6 0 3\nS7 S7 0 4" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2\nS6 S6 0 4\nS7 S7 0 NA" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2\nS6 S6 0 NA\nS7 S7 0 3" ) ;

		// Run the tests.
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "test #" << i << "\n" ;
			std::istringstream istr( test_data[i] ) ;
			genfile::TraditionalStrictCohortIndividualSource source( istr ) ;
			genfile::CategoricalCrossCohortCovariateValueMapping mapping( "column1" ) ;
			mapping.add_source( source ) ;
			std::size_t number_of_missing_values = std::count( test_data[i].begin(), test_data[i].end(), 'N' ) ;
			TEST_ASSERT( mapping.get_number_of_missing_values() == number_of_missing_values ) ;
			TEST_ASSERT( mapping.get_number_of_distinct_mapped_values() == 4 ) ;
			TEST_ASSERT( mapping.get_unmapped_value(1) == genfile::CohortIndividualSource::Entry( 1 ) ) ;
			TEST_ASSERT( mapping.get_unmapped_value(2) == genfile::CohortIndividualSource::Entry( 2 ) ) ;
			TEST_ASSERT( mapping.get_unmapped_value(3) == genfile::CohortIndividualSource::Entry( 3 ) ) ;
			TEST_ASSERT( mapping.get_unmapped_value(4) == genfile::CohortIndividualSource::Entry( 4 ) ) ;
		}
	}

	void multiple_cohort_test() {
		std::cerr << "CategoricalValueMappingTester::multiple_cohort_test()\n" ;
		
		std::vector< std::vector< std::string > > test_data ;
		std::vector< std::vector< genfile::CohortIndividualSource::Entry > > expected ;
	
		// 1. test that a few cohorts with the same value only give one level
		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1 ) ;

		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1 ) ;

		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1 ) ;

		// 1. test a few cohorts with differing values
		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 2" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1 ) ;
		expected.back().push_back( 2 ) ;

		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 2" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 2" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1 ) ;
		expected.back().push_back( 2 ) ;

		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 2" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 D\nS1 S1 0 3" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1 ) ;
		expected.back().push_back( 2 ) ;
		expected.back().push_back( 3 ) ;

		// Run the tests.
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "test #" << i << "\n" ;
			genfile::CategoricalCrossCohortCovariateValueMapping mapping( "column1" ) ;

			for( std::size_t j = 0; j < test_data[i].size(); ++j ) {
				std::istringstream istr( test_data[i][j] ) ;
				genfile::TraditionalStrictCohortIndividualSource source( istr ) ;
				mapping.add_source( source ) ;
			}
			TEST_ASSERT( mapping.get_number_of_distinct_mapped_values() == expected[i].size() ) ;
			for( std::size_t j = 0; j < expected[i].size(); ++j ) {
				TEST_ASSERT( mapping.get_unmapped_value( int(j+1) ) == expected[i][j] ) ;
			}
		}

	}
} ;

struct ContinuousValueMappingTester
{
	ContinuousValueMappingTester() {
		number_of_mapped_values_test() ;
		level_value_test() ;
		multiple_cohort_test() ;
	}
	
	// Test that with a single source we get the correct number of levels.
	void number_of_mapped_values_test() {
		std::cerr << "ContinuousValueMappingTester::number_of_mapped_values_test()\n" ;
		std::vector< std::string > test_data ;
		std::vector< std::size_t > expected_number_of_levels ;
		
		// 0 -entry tests
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C" ) ;
		expected_number_of_levels.push_back( 0 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 P" ) ;
		expected_number_of_levels.push_back( 0 ) ;

		// Continuous covariate tests
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 0.1" ) ;
		expected_number_of_levels.push_back( 1 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 0.1\nS2 S2 0 0.2" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 0.1\nS2 S2 0 0.3" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 0.1\nS2 S2 0 0.2\nS3 S3 0 0.3" ) ;
		expected_number_of_levels.push_back( 3 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 0.1\nS2 S2 0 0.2\nS3 S3 0 0.2" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 0.1\nS2 S2 0 0.2\nS3 S3 0 0.2\n0 0 0 0.1" ) ;
		expected_number_of_levels.push_back( 2 ) ;

		// Continuous phenotype tests
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 P\nS1 S1 0 0.1" ) ;
		expected_number_of_levels.push_back( 1 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 P\nS1 S1 0 0.1\nS2 S2 0 0.2" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 P\nS1 S1 0 0.1\nS2 S2 0 0.3" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 P\nS1 S1 0 0.1\nS2 S2 0 0.2\nS3 S3 0 0.3" ) ;
		expected_number_of_levels.push_back( 3 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 P\nS1 S1 0 0.1\nS2 S2 0 0.2\nS3 S3 0 0.2" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 P\nS1 S1 0 0.1\nS2 S2 0 0.2\nS3 S3 0 0.2\nS4 S4 0 0.1" ) ;
		expected_number_of_levels.push_back( 2 ) ;
		
		// For each of the above, try again with a number of missing values added.
		{
			std::size_t const N = test_data.size() ;
			for( std::size_t i = 0; i < N; ++i ) {
				std::string data = test_data[i] ;
				for( std::size_t j = 1; j < 4; ++j ) {
					data += "\nSM" + genfile::string_utils::to_string( j ) + " SM" + genfile::string_utils::to_string( j ) + " 0 NA" ;
				}
				test_data.push_back( data ) ;	
				expected_number_of_levels.push_back( expected_number_of_levels[i] ) ;
			}
		}
		
		// Run the tests.
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "test #" << i << "\n" ;
			std::istringstream istr( test_data[i] ) ;
			genfile::TraditionalStrictCohortIndividualSource source( istr ) ;
			genfile::ContinuousVariableCrossCohortCovariateValueMapping mapping( "column1" ) ;
			mapping.add_source( source ) ;
			TEST_ASSERT( mapping.get_number_of_distinct_mapped_values() == expected_number_of_levels[i] ) ;
		}
	}

	void level_value_test() {
		std::cerr << "ContinuousValueMappingTester::level_value_test()\n" ;
		
		std::vector< std::string > test_data ;
		std::vector< genfile::CohortIndividualSource::Entry > expecteds ;
		
		// In all cases the levels should be NA,1,2,3,4.
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 3\nS4 S4 0 4\nS5 S5 0 NA" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 NA\nS2 S2 0 4\nS3 S3 0 3\nS4 S4 0 2\nS5 S5 0 1" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 4\nS2 S2 0 1\nS3 S3 0 2\nS4 S4 0 3\nS5 S5 0 NA" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2" ) ;

		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 3\nS4 S4 0 4\nS5 S5 0 NA\nS6 S6 0 1" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 NA\nS2 S2 0 4\nS3 S3 0 3\nS4 S4 0 2\nS5 S5 0 1\nS6 S6 0 2" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 4\nS2 S2 0 1\nS3 S3 0 2\nS4 S4 0 3\nS5 S5 0 NA\nS6 S6 0 3" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2\nS6 S6 0 4" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2\nS6 S6 0 NA" ) ;

		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 3\nS4 S4 0 4\nS5 S5 0 NA\nS6 S6 0 1\nS7 S7 0 2" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 NA\nS2 S2 0 4\nS3 S3 0 3\nS4 S4 0 2\nS5 S5 0 1\nS6 S6 0 2\nS7 S7 0 1" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 4\nS2 S2 0 1\nS3 S3 0 2\nS4 S4 0 3\nS5 S5 0 NA\nS6 S6 0 3\nS7 S7 0 4" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2\nS6 S6 0 4\nS7 S7 0 NA" ) ;
		test_data.push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 3\nS2 S2 0 NA\nS3 S3 0 4\nS4 S4 0 1\nS5 S5 0 2\nS6 S6 0 NA\nS7 S7 0 3" ) ;

		// Run the tests.
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "test #" << i << "\n" ;
			std::istringstream istr( test_data[i] ) ;
			genfile::TraditionalStrictCohortIndividualSource source( istr ) ;
			genfile::ContinuousVariableCrossCohortCovariateValueMapping mapping( "column1" ) ;
			mapping.add_source( source ) ;
			std::size_t number_of_missing_values = std::count( test_data[i].begin(), test_data[i].end(), 'N' ) ;
			TEST_ASSERT( mapping.get_number_of_missing_values() == number_of_missing_values ) ;
			TEST_ASSERT( mapping.get_unmapped_value( 1.0 ) == genfile::CohortIndividualSource::Entry( 1.0 ) ) ;
			TEST_ASSERT( mapping.get_unmapped_value( 2.0 ) == genfile::CohortIndividualSource::Entry( 2.0 ) ) ;
			TEST_ASSERT( mapping.get_unmapped_value( 3.0 ) == genfile::CohortIndividualSource::Entry( 3.0 ) ) ;
			TEST_ASSERT( mapping.get_unmapped_value( 4.0 ) == genfile::CohortIndividualSource::Entry( 4.0 ) ) ;
		}
	}

	void multiple_cohort_test() {
		std::cerr << "ContinuousValueMappingTester::multiple_cohort_test()\n" ;
		
		std::vector< std::vector< std::string > > test_data ;
		std::vector< std::vector< genfile::CohortIndividualSource::Entry > > expected ;
	
		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1.0 ) ;

		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 2" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1.0 ) ;
		expected.back().push_back( 2.0 ) ;

		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 2" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 3" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 1.0 ) ;
		expected.back().push_back( 2.0 ) ;
		expected.back().push_back( 3.0 ) ;

		// Run the tests.
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "test #" << i << "\n" ;
			genfile::ContinuousVariableCrossCohortCovariateValueMapping mapping( "column1" ) ;

			for( std::size_t j = 0; j < test_data[i].size(); ++j ) {
				std::istringstream istr( test_data[i][j] ) ;
				genfile::TraditionalStrictCohortIndividualSource source( istr ) ;
				mapping.add_source( source ) ;
			}
			TEST_ASSERT( mapping.get_number_of_distinct_mapped_values() == expected[i].size() ) ;
			for( std::size_t j = 0; j < expected[i].size(); ++j ) {
				TEST_ASSERT( mapping.get_unmapped_value( double(j+1) ) == genfile::CohortIndividualSource::Entry( expected[i][j] )) ;
			}
		}

	}
} ;



struct NormalisingValueMappingTester
{
	NormalisingValueMappingTester() {
		multiple_cohort_test() ;
	}
	
	void multiple_cohort_test() {
		std::cerr << "NormalisingValueMappingTester::multiple_cohort_test()\n" ;
		
		std::vector< std::vector< std::string > > test_data ;
		std::vector< std::vector< genfile::CohortIndividualSource::Entry > > expected ;
	
		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1.0" ) ;
		expected.resize( expected.size() + 1 ) ;
		expected.back().push_back( 0.0 ) ;
	
		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1.0" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 2.0" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 3.0" ) ;
		expected.resize( expected.size() + 1 ) ;
		// variance and sd are 1.
		expected.back().push_back( -1.0 ) ;
		expected.back().push_back( 0.0 ) ;
		expected.back().push_back( 1.0 ) ;

		test_data.resize( test_data.size() + 1 ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 1.0" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 2.0" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 3.0" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 4.0" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 5.0" ) ;
		test_data.back().push_back( "ID_1 ID_2 missing column1\n0 0 0 C\nS1 S1 0 6.0" ) ;
		expected.resize( expected.size() + 1 ) ;
		{
			// variance is 3.5
			double sd = std::sqrt( 3.5 ) ;
			expected.back().push_back( -2.5 / sd ) ;
			expected.back().push_back( -1.5 / sd ) ;
			expected.back().push_back( -0.5 / sd ) ;
			expected.back().push_back( 0.5 / sd ) ;
			expected.back().push_back( 1.5 / sd ) ;
			expected.back().push_back( 2.5 / sd ) ;
		}
		// Run the tests.
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "test #" << i << "\n" ;
			genfile::NormalisingCrossCohortCovariateValueMapping mapping( "column1" ) ;

			for( std::size_t j = 0; j < test_data[i].size(); ++j ) {
				std::istringstream istr( test_data[i][j] ) ;
				genfile::TraditionalStrictCohortIndividualSource source( istr ) ;
				mapping.add_source( source ) ;
			}
			TEST_ASSERT( mapping.get_number_of_distinct_mapped_values() == expected[i].size() ) ;
			for( std::size_t j = 0; j < expected[i].size(); ++j ) {
				std::cerr << mapping.get_mapped_value( int(j+1) ) << " " << genfile::CohortIndividualSource::Entry( expected[i][j] ) << ".\n" ;
				TEST_ASSERT( std::abs( mapping.get_mapped_value( int(j+1) ).as< double >() - genfile::CohortIndividualSource::Entry( expected[i][j] ).as< double >() ) < 0.00000001 ) ;
			}
		}
	}
} ;
