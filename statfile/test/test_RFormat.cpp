#include <iostream>
#include <utility>
#include <map>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "test_case.hpp"
#include "statfile/RFormatStatSink.hpp"
#include "statfile/RFormatStatSource.hpp"

namespace globals {
	std::string data =
	"# This is a file\n"
	"# These lines are comments, which should be ignored.\n"
	"#\n"
	"index Column1 Column2 Column3 0 1 2 3 4 5 6 7 8 9\n"
	"0 Hello hello hello 0 0 0 0 0 0 0 0 0 0\n" 
	"1 H e L 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n" 
	"100 Here's a test 100.0 100000.254 576.22 -0.1 -100.001 100.001 1000000 4294967296 9007199254740992 -9007199254740992\n"
	"110 A B C 0.00001 9999.99999 0.00001 9999.99999 0.00001 9999.99999 0.00001 9999.99999 0.00001 9999.99999\n" ;

	int number_of_data_columns = 10 ;
	int number_of_rows = 4 ;
}

void copy_data_to_file( std::string const& data, std::string const& filename ) {
	std::ofstream file( filename.c_str() ) ;
	TEST_ASSERT( file.is_open() ) ;
	file << data ;
	file.close() ;
}

AUTO_TEST_CASE( test_RFormat ) 
{
	std::string filename = std::tmpnam(0) ;
	std::string filename2 = std::tmpnam(0) ;

	std::cout << "First filename: \"" << filename << "\".  Second filename: \"" << filename2 << "\".\n" ;

	copy_data_to_file( globals::data, filename ) ;

	{
		statfile::RFormatStatSource stat_source_1( filename ) ;
		statfile::RFormatStatSink stat_sink_1( filename2 ) ;
		stat_sink_1.set_descriptive_text( stat_source_1.get_descriptive_text()) ;
		int index ;
		std::string column1, column2, column3 ;
		std::vector< double > doubles( 10 ) ;

		stat_sink_1.add_columns( stat_source_1.column_names() ) ;

		while( stat_source_1 >> index ) {
			stat_source_1 >> column1 >> column2 >> column3 ;
			for( std::size_t i = 0; i < doubles.size(); ++i ) {
				stat_source_1 >> doubles[i] ;
			}

			stat_sink_1 << index << column1 << column2 << column3 ;
			for( std::size_t i = 0; i < doubles.size(); ++i ) {
				stat_sink_1 << doubles[i] ;
			}
		}
	}
	
	{
		statfile::RFormatStatSource stat_source_1( filename ) ;
		statfile::RFormatStatSource stat_source_2( filename2 ) ;
		TEST_ASSERT( stat_source_1.get_descriptive_text() == stat_source_2.get_descriptive_text() ) ;

		std::vector< double > doubles1( 10 ) ;
		std::vector< double > doubles2( 10 ) ;
	
		int index1 ;
		std::string column11, column12, column13 ;
		int index2 ;
		std::string column21, column22, column23 ;
	
		double epsilon = 0.00000000001 ;
	
		while(( stat_source_1 >> index1 ) && ( stat_source_2 >> index2 )) {
			stat_source_1 >> column11 >> column12 >> column13 ;
			stat_source_2 >> column21 >> column22 >> column23 ;
			for( std::size_t i = 0; i < 10; ++i ) {
				stat_source_1 >> doubles1[i] ;
				stat_source_2 >> doubles2[i] ;
			}
			TEST_ASSERT( stat_source_1 ) ;
			TEST_ASSERT( stat_source_2 ) ;
			std::cout << "Comparing lines with index " << index1 << ", " << index2 << "...\n" ;
			TEST_ASSERT( index1 == index2 ) ;
			TEST_ASSERT( column11 == column21 ) ;
			TEST_ASSERT( column12 == column22 ) ;
			TEST_ASSERT( column13 == column23 ) ;
			TEST_ASSERT( doubles1.size() == doubles2.size() ) ;
			for( std::size_t i = 0; i < doubles1.size(); ++i ) {
				TEST_ASSERT( std::abs( doubles1[i] - doubles2[i]) < epsilon ) ;
			}
		}
	}
	
	std::cout << "Success.\n" ;
}

AUTO_TEST_MAIN
{
	test_RFormat() ;
}
