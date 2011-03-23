#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "test_case.hpp"
#include "genfile/TraditionalStrictCohortIndividualSource.hpp"
#include "genfile/CategoricalCohortIndividualSource.hpp"
#include <cassert>

typedef boost::function< genfile::CohortIndividualSource::UniquePtr ( std::istream&, std::string const& ) > SourceConstructor ;
genfile::CohortIndividualSource::UniquePtr create_strict_source( std::istream& stream, std::string const& missing_values ) {
	return genfile::CohortIndividualSource::UniquePtr( new genfile::TraditionalStrictCohortIndividualSource( stream, missing_values )) ;
}

genfile::CohortIndividualSource::UniquePtr create_categorical_source( std::istream& stream, std::string const& missing_values ) {
	return genfile::CohortIndividualSource::UniquePtr( new genfile::CategoricalCohortIndividualSource( stream, missing_values )) ;
}

struct TraditionalStrictCohortIndividualSourceTester
{
	TraditionalStrictCohortIndividualSourceTester() {
		test_column_headers() ;
		test_mandatory_columns() ;
		test_continuous_columns() ;
		test_discrete_columns() ;
		test_binary_columns() ;
		test_whitespace() ;
		test_sample_ids() ;
		test_missing_values() ;
	}
	
	void test_column_headers() {
		std::cerr << "test_column_headers\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
	
		test_data.push_back( "ID_1 ID_2 missing" ) ;
		expected_success.push_back( false ) ; // no column types row
		test_data.push_back( "ID_1 ID_2 missing\n" ) ;
		expected_success.push_back( false ) ; // no column types row
		test_data.push_back( "ID_1 ID_2 missing\n0" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing \n0 0 0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing \n0 0 0\nSAMPLE1 SAMPLE2 0.0" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing \n0 0 0\nSAMPLE1 SAMPLE2 0.0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing \n0 0 0\nSAMPLE1 ASAMPLE1 0.0\nSAMPLE2 ASAMPLE2 10.0" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing \n0 0 0\nSAMPLE1 ASAMPLE1 0.0\nSAMPLE2 ASAMPLE2 10.0\n" ) ;
		expected_success.push_back( true ) ;
		
		run_test( test_data, expected_success, &create_strict_source ) ;		
		run_test( test_data, expected_success, &create_categorical_source ) ;		
	}

	void test_mandatory_columns() {
		std::cerr << "test_mandatory_columns\n" ;
		
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\na a a\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\na a 1.0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\na a 25.0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nhello there 25.0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nhello there number\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nWT_1_556_kj WXXX&&&^@@ 25.0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nWT_1_556_kj WXXX&&&^@@ &&&6\n" ) ;
		expected_success.push_back( false ) ;

		run_test( test_data, expected_success, &create_strict_source ) ;		
		run_test( test_data, expected_success, &create_categorical_source ) ;		
	}

	void test_continuous_columns() {
		std::cerr << "test_continuous_columns\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
		// continuous covariates
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 C\nS1 S1 0.0 1" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 C\nS1 S1 0.0 1.0" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 C\nS1 S1 0.0 1234.5" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 C\nS1 S1 0.0 hello" ) ;
		expected_success.push_back( false ) ;
		
		// continuous phenotype
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 P\nS1 S1 0.0 1" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 P\nS1 S1 0.0 1.0" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 P\nS1 S1 0.0 1234.5" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 P\nS1 S1 0.0 hello" ) ;
		expected_success.push_back( false ) ;

		run_continuous_value_test( test_data, expected_success, "acolumn", &create_strict_source ) ;
		run_continuous_value_test( test_data, expected_success, "acolumn", &create_categorical_source ) ;
	}
	
	void test_discrete_columns() {
		std::cerr << "test_discrete_columns\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success_strict ;

		// positive integral values are fine
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 D\nS1 S1 0.0 1" ) ;
		expected_success_strict.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 D\nS1 S1 0.0 1234" ) ;
		expected_success_strict.push_back( true ) ;

		// Non-integral values should fail
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 D\nS1 S1 0.0 1.5" ) ;
		expected_success_strict.push_back( false ) ;

		// Non-numerical values should fail
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 D\nS1 S1 0.0 hello" ) ;
		expected_success_strict.push_back( false ) ;

		// zero or -ve values should fail
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 D\nS1 S1 0.0 0" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 D\nS1 S1 0.0 -1" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 D\nS1 S1 0.0 -1234" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 D\nS1 S1 0.0 -1.1" ) ;
		expected_success_strict.push_back( false ) ;

		// All of the above should be ok for categorical source.
		std::vector< bool > expected_success_categorical( expected_success_strict.size(), true ) ;
		
		run_discrete_value_test( test_data, expected_success_strict, "acolumn", &create_strict_source ) ;
		run_categorical_value_test( test_data, expected_success_categorical, "acolumn", &create_categorical_source ) ;
	}

	void test_binary_columns() {
		std::cerr << "test_binary_columns\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success_strict ;
		std::vector< bool > expected_success_categorical ;
		// binary phenotype.  Value should be zero or one.
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 0" ) ;
		expected_success_strict.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 1" ) ;
		expected_success_strict.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 -1" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 1.1" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 2" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 10" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 1.0" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 1234.5" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 h" ) ;
		expected_success_strict.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 hello" ) ;
		expected_success_strict.push_back( false ) ;

		// up to now, strict and categorical sources should behave the same.
		expected_success_categorical = expected_success_strict ;

		// We now explore the differences.
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 case" ) ;
		expected_success_strict.push_back( false ) ;
		expected_success_categorical.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing acolumn\n0 0 0 B\nS1 S1 0.0 control" ) ;
		expected_success_strict.push_back( false ) ;
		expected_success_categorical.push_back( true ) ;

		run_discrete_value_test( test_data, expected_success_strict, "acolumn", &create_strict_source ) ;
		run_discrete_value_test( test_data, expected_success_categorical, "acolumn", &create_categorical_source ) ;
	}

	void test_number_of_columns_requirement() {
		std::cerr << "test_number_of_columns_requirement\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
	
		test_data.push_back( "ID_1 ID_2 missing" ) ;
		expected_success.push_back( false ) ; // no column types row
		test_data.push_back( "ID_1 ID_2 missing\n" ) ;
		expected_success.push_back( false ) ; // no column types row

		test_data.push_back( "ID_1 ID_2 missing\n0 0 0" ) ;
		expected_success.push_back( true ) ;

		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\n" ) ;
		expected_success.push_back( true ) ;

		test_data.push_back( "ID_1 ID_2 missing \n0 0 0\n" ) ;
		expected_success.push_back( true ) ;
		
		run_test( test_data, expected_success, &create_strict_source ) ;		
		run_test( test_data, expected_success, &create_categorical_source ) ;		
	}

	void test_several_rows() {
		std::cerr << "test_several_rows\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
	
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0 12" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0 12\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0 12\nS2 S2 130 11.3 1 11.56" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0 12\nS2 S2 130 11.3 1 11.56\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0\nS2 S2 130 11.3 1 11.56\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0\nS2 S2 130 11.3 1" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0\nS2 S2 130 11.3 1\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0 12\nS2 S2 130 11.3 1 11.56 13" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 12 12 0 12\nS2 S2 130 11.3 1 11.56 13\n" ) ;
		expected_success.push_back( false ) ;
		
		run_test( test_data, expected_success, &create_strict_source ) ;		
		run_test( test_data, expected_success, &create_categorical_source ) ;		
	}
	
	void test_whitespace() {
		std::cerr << "test_whitespace\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
	
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 0.0 12 12 0 12\nS2 S2 0.5 130 11.3 1 11.56" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 0.0 12 12 0 12\nS2 S2 0.5 130 11.3 1 11.56\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 0.0 12 12 0 12\nS2 S2 0.5 130 11.3 1 11.56 " ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 0.0 12 12 0 12\nS2 S2 0.5 130 11.3 1 11.56 \n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 0.0 12 12 0 12\nS2 S2 0.5 130 11.3 1 11.56  " ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 0.0 12 12 0 12\nS2 S2 0.5 130 11.3 1 11.56  \n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1\tID_2\tmissing\ta\tb\tc\td\n0\t0\t0\tD\tC\tB\tP\nS1 S1 0.0 12\t12\t0\t12\nS2 S2 0.5 130\t11.3\t1\t11.56" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1\tID_2\tmissing\ta\tb\tc\td\n0\t0\t0\tD\tC\tB\tP\nS1 S1 0.0 12\t12\t0\t12\nS2 S2 0.5 130\t11.3\t1\t11.56\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1\tID_2\tmissing\ta\tb\tc\td\n0\t0\t0\tD\tC\tB\tP\nS1 S1 0.0 12\t12\t0\t12\nS2 S2 0.5 130\t11.3\t1\t11.56\t" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1\tID_2\tmissing\ta\tb\tc\td\n0\t0\t0\tD\tC\tB\tP\nS1 S1 0.0 12\t12\t0\t12\nS2 S2 0.5 130\t11.3\t1\t11.56\t\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1\tID_2\tmissing\ta\tb\tc\td\n0\t0\t0\tD\tC\tB\tP\nS1 S1 0.0 12\t12\t0\t12\nS2 S2 0.5 130\t11.3\t1\t11.56\t\t" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1\tID_2\tmissing\ta\tb\tc\td\n0\t0\t0\tD\tC\tB\tP\nS1 S1 0.0 12\t12\t0\t12\nS2 S2 0.5 130\t11.3\t1\t11.56\t\t\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1  ID_2  missing  a   b    c  d\n 0  0  0    D  C  B  P   \n    \tS1\tS1  0.0   \t 12  \t 12  0  \t12\nS2 S2 0.5 130  11.3   1  11.56" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1  ID_2  missing  a   b    c  d\n 0  0  0    D  C  B  P   \n    \tS1\tS1  \t0.0   \t 12  \t 12  0  \t12\nS2 S2 0.5 130  11.3   1  11.56 \n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1  ID_2  missing  a   b    c  d\n 0  0  0    D  C  B  P   \n    \tS1\tS1   \t0.0\t  \t 12  \t 12  0  \t12\nS2 S2 0.5 130  11.3   1  11.56\t\n" ) ;
		expected_success.push_back( true ) ;

		test_data.push_back( "ID_1  ID_2  missing  a   b    c  d\n 0  0 \n0    D  C  B  P   \n    \t  S1\tS1  0.0   \t 12  \t 12  0  \t12\nS2 S2 0.5 130  11.3   1  11.56\t\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1  \nID_2  missing  a   b    c  d\n 0  0  0    D  C  B  P   \n    \t S1\tS1  \t0.0   \t 12  \t 12  0  \t12\nS2 S2 0.5 130  11.3   1  11.56\t\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a\nb c d\n0 0 0 D C B P\nS1 S1 0.0 12 12 0 12\nS2 S2 0.5 130 11.3 1 11.56" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a b c d\n0 0 0 D C B P\nS1 S1 0.0 12 12 0 12\nS2 S2 0.5 130 11.3 1 11.56\n\n" ) ;
		expected_success.push_back( false ) ;

		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "case " << i << "..." ;
			try {
				{
					std::cerr << "strict..." ;
					std::istringstream is( test_data[i] ) ;
					genfile::TraditionalStrictCohortIndividualSource source( is ) ;
					TEST_ASSERT( expected_success[i] ) ;

					TEST_ASSERT( source.get_entry( 0, "a" ).as< int >() == 12 ) ;
					TEST_ASSERT( source.get_entry( 0, "b" ).as< double >() == 12 ) ;
					TEST_ASSERT( source.get_entry( 0, "c" ).as< int >() == 0 ) ;
					TEST_ASSERT( source.get_entry( 0, "d" ).as< double >() == 12 ) ;
					TEST_ASSERT( source.get_entry( 1, "a" ).as< int >() == 130 ) ;
					TEST_ASSERT( source.get_entry( 1, "b" ).as< double >() == 11.3 ) ;
					TEST_ASSERT( source.get_entry( 1, "c" ).as< int >() == 1 ) ;
					TEST_ASSERT( source.get_entry( 1, "d" ).as< double >() == 11.56 ) ;
				}
				{
					std::cerr << "categorical..." ;
					std::istringstream is( test_data[i] ) ;
					genfile::CategoricalCohortIndividualSource source( is ) ;
					TEST_ASSERT( expected_success[i] ) ;

					TEST_ASSERT( source.get_entry( 0, "a" ).as< std::string >() == "12" ) ;
					TEST_ASSERT( source.get_entry( 0, "b" ).as< double >() == 12 ) ;
					TEST_ASSERT( source.get_entry( 0, "c" ).as< int >() == 0 ) ;
					TEST_ASSERT( source.get_entry( 0, "d" ).as< double >() == 12 ) ;
					TEST_ASSERT( source.get_entry( 1, "a" ).as< std::string >() == "130" ) ;
					TEST_ASSERT( source.get_entry( 1, "b" ).as< double >() == 11.3 ) ;
					TEST_ASSERT( source.get_entry( 1, "c" ).as< int >() == 1 ) ;
					TEST_ASSERT( source.get_entry( 1, "d" ).as< double >() == 11.56 ) ;
				}
			}
			catch( genfile::InputError const& e ) {
				TEST_ASSERT( !expected_success[i] ) ;
			}
		}
		std::cerr << "done\n" ;

	}
	
	void test_sample_ids() {
		std::cerr << "test_sample_ids\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
	
		// We should not be allowed repeated entries in ID_1, regardless of ID_2
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nS1 S1 0.0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nS1 S1A 0.0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nS1 S1A 0.0\nS1 S1A 0.0\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nS1 S1A 0.0\nS2 S1A 0.0\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nS1 S1A 0.0\nS1 S2A 0.0\n" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nS1 S1A 0.0\nS2 S2A 0.0\n" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nS1 S1A 0.0\nS2 S2A 0.0\nS1 S1 0.0" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing\n0 0 0\nS1 S1A 0.0\nS2 S2A 0.0\nS1 S1A 0.0" ) ;
		expected_success.push_back( false ) ;

		run_test( test_data, expected_success, &create_strict_source ) ;		
		run_test( test_data, expected_success, &create_categorical_source ) ;		
	}
	
	void test_missing_values() {
		std::vector< std::string > missing_values ;
		
		// Single missing values
		missing_values.push_back( "NA" ) ;
		missing_values.push_back( " NA " ) ;
		missing_values.push_back( "NA2" ) ;
		missing_values.push_back( "2NA" ) ;
		missing_values.push_back( "-999" ) ;
		missing_values.push_back( "999" ) ;
		missing_values.push_back( "OhMyGoodness" ) ;
		missing_values.push_back( "-&8&@*" ) ;
		
		// multiple missing values
		missing_values.push_back( "NA,NA2" ) ;
		missing_values.push_back( "NA, NA2" ) ;
		missing_values.push_back( " NA,    NA2  " ) ;
		missing_values.push_back( "NA , NA2" ) ;
		missing_values.push_back( "NA, NA2, NA3" ) ;
		missing_values.push_back( "NA, NA2, NA3, NA4" ) ;
		missing_values.push_back( " NA  , NA2, NA3 , NA4 " ) ;

		for( std::size_t i = 0; i < missing_values.size(); ++i ) {
			// get a list of actual missing values.
			std::size_t pos1 = missing_values[i].find_first_not_of( " ," ) ;
			std::size_t pos2 = missing_values[i].find_first_of( " ,", pos1 ) ;
			std::vector< std::string > missing_value_list ;
			while( pos1 != std::string::npos ) {
				missing_value_list.push_back( missing_values[i].substr( pos1, pos2 - pos1 )) ;
				pos1 = ( pos2 == std::string::npos ) ? std::string::npos : missing_values[i].find_first_not_of( " ,", pos2 ) ;
				pos2 = ( pos1 == std::string::npos ) ? std::string::npos : missing_values[i].find_first_of( " ,", pos1 ) ;
			}
			for( std::size_t j = 0; j < missing_value_list.size(); ++j ) {
				test_missing_value_parsing( missing_values[i], missing_value_list[j] ) ;
				test_missing_value_retrieval( missing_values[i], missing_value_list[j] ) ;
			}
		}
	}

	void test_missing_value_parsing( std::string const& missing_values, std::string const& missing_value ) {
		std::cerr << "test_missing_value_parsing( \"" << missing_values << "\", \"" << missing_value << "\" )\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
	
		// Test that we get missing values in missing columns
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 D\nS1 S1 0.0 " + missing_value ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 D\nS1 S1 0.0 " + missing_value + "\nS2 S2 0.0 1" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 B\nS1 S1 0.0 " + missing_value + "\nS2 S2 0.0 1" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 B\nS1 S1 0.0 " + missing_value ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 C\nS1 S1 0.0 " + missing_value ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 C\nS1 S1 0.0 " + missing_value + "\nS2 S2 0.0 10.0" ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 P\nS1 S1 0.0 " + missing_value ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 P\nS1 S1 0.0 " + missing_value + "\nS2 S2 0.0 10.0" ) ;
		expected_success.push_back( true ) ;

		// Check that the file parses
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			try {
				{
					std::istringstream is( test_data[i] ) ;
					genfile::TraditionalStrictCohortIndividualSource source( is, missing_values ) ;
					TEST_ASSERT( expected_success[i] ) ;
				}
				{
					std::istringstream is( test_data[i] ) ;
					genfile::CategoricalCohortIndividualSource source( is, missing_values ) ;
					TEST_ASSERT( expected_success[i] ) ;
				}
			}
			catch( genfile::InputError const& e ) {
				TEST_ASSERT( !expected_success[i] ) ;
			}
		}
	}

	void test_missing_value_retrieval( std::string const& missing_values, std::string const& missing_value ) {
		std::cerr << "test_missing_value_retrieval( \"" << missing_values << "\", \"" << missing_value << "\" )\n" ;
		std::vector< std::string > test_data ;
		std::vector< bool > expected_success ;
	
		// Test that we get missing values put in missing columns
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 D\nS1 S1 0.0 " + missing_value ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 D\nS1 S1 0.0 1" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 B\nS1 S1 0.0 " + missing_value ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 B\nS1 S1 0.0 1" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 C\nS1 S1 0.0 " + missing_value ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 C\nS1 S1 0.0 10.0" ) ;
		expected_success.push_back( false ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 P\nS1 S1 0.0 " + missing_value ) ;
		expected_success.push_back( true ) ;
		test_data.push_back( "ID_1 ID_2 missing a_column\n0 0 0 P\nS1 S1 0.0 10.0" ) ;
		expected_success.push_back( false ) ;

		// Check that the value is retrived as missing
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			{
				std::istringstream is( test_data[i] ) ;
				genfile::TraditionalStrictCohortIndividualSource source( is, missing_value ) ;
				if( source.get_entry( 0, "a_column" ).is_missing() ) {
					TEST_ASSERT( expected_success[i] ) ;
				}
				else {
					TEST_ASSERT( !expected_success[i] ) ;
				}
			}
			{
				std::istringstream is( test_data[i] ) ;
				genfile::CategoricalCohortIndividualSource source( is, missing_value ) ;
				if( source.get_entry( 0, "a_column" ).is_missing() ) {
					TEST_ASSERT( expected_success[i] ) ;
				}
				else {
					TEST_ASSERT( !expected_success[i] ) ;
				}
			}
		}
		
	}

	void run_test( std::vector< std::string > const& test_data, std::vector< bool > const& expected_success, SourceConstructor source_constructor ) const {
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "case " << i << "..." ;
			try {
				std::istringstream is( test_data[i] ) ;
				genfile::CohortIndividualSource::UniquePtr source( source_constructor( is, "NA" )) ;
				TEST_ASSERT( expected_success[i] ) ;
			}
			catch( genfile::InputError const& e ) {
				TEST_ASSERT( !expected_success[i] ) ;
			}
		}
		std::cerr << "done\n" ;
	}

	void run_discrete_value_test(
		std::vector< std::string > const& test_data,
		std::vector< bool > const& expected_success,
		std::string const& column_name,
		SourceConstructor source_constructor
	) const {
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "case " << i << "..." ;
			try {
				std::istringstream is( test_data[i] ) ;
				genfile::CohortIndividualSource::UniquePtr source = source_constructor( is, "NA" ) ;
				source->get_entry( 0, column_name ).as< int >() ;
				source->get_entry( 0, column_name ).as< double >() ;
				TEST_ASSERT( expected_success[i] ) ;
			}
			catch( genfile::InputError const& e ) {
				TEST_ASSERT( !expected_success[i] ) ;
			}
			catch( boost::bad_get const& e ) {
				TEST_ASSERT( !expected_success[i] ) ;
			}
		}
		std::cerr << "done\n" ;
	}

	void run_categorical_value_test(
		std::vector< std::string > const& test_data,
		std::vector< bool > const& expected_success,
		std::string const& column_name,
		SourceConstructor source_constructor
	) const {
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "case " << i << "..." ;
			try {
				std::istringstream is( test_data[i] ) ;
				genfile::CohortIndividualSource::UniquePtr source = source_constructor( is, "NA" ) ;
				source->get_entry( 0, column_name ).as< std::string >() ;
				TEST_ASSERT( expected_success[i] ) ;
			}
			catch( genfile::InputError const& e ) {
				TEST_ASSERT( !expected_success[i] ) ;
			}
			catch( boost::bad_get const& e ) {
				TEST_ASSERT( !expected_success[i] ) ;
			}
		}
		std::cerr << "done\n" ;
	}

	void run_continuous_value_test(
		std::vector< std::string > const& test_data,
		std::vector< bool > const& expected_success,
		std::string const& column_name,
		SourceConstructor source_constructor
	) const {
		for( std::size_t i = 0; i < test_data.size(); ++i ) {
			std::cerr << "case " << i << "..." ;
			try {
				std::istringstream is( test_data[i] ) ;
				genfile::CohortIndividualSource::UniquePtr source = source_constructor( is, "NA" ) ;
				source->get_entry( 0, column_name ).as< double >() ;
				TEST_ASSERT( expected_success[i] ) ;
			}
			catch( genfile::InputError const& e ) {
				TEST_ASSERT( !expected_success[i] ) ;
			}
		}
		std::cerr << "done\n" ;
	}
} ;

AUTO_TEST_CASE( all_tests ) {
	TraditionalStrictCohortIndividualSourceTester tester ;
}

AUTO_TEST_MAIN {
	all_tests() ;
}
