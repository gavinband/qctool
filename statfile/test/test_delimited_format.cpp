
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <utility>
#include <map>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "test_case.hpp"
#include "statfile/DelimitedStatSource.hpp"

namespace globals {
	namespace {
		std::string data =
		"# This is a file\n"
		"# These lines are comments, which should be ignored.\n"
		"#\n"
		"index\tColumn1\tColumn2\tColumn3\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\n"
		"0\tHello\thello\thello\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
		"1\tH\te\tL\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n"
		"2\tHere's\ta\ttest\t100.0\t100000.254\t576.22\t-0.1\t-100.001\t100.001\t1000000\t4294967296\t9007199254740992\t-9007199254740992\n"
		"3\tA\tB\tC\t0.00001\t9999.99999\t0.00001\t9999.99999\t0.00001\t9999.99999\t0.00001\t9999.99999\t0.00001\t9999.99999\n" ;

		std::string quoted_data =
		"# This is a file\n"
		"# These lines are comments, which should be ignored.\n"
		"#\n"
		"\"index\"\t\"Column1\"\t\"Column2\"\t\"Column3\"\t\"0\"\t\"1\"\t\"2\"\t\"3\"\t\"4\"\t\"5\"\t\"6\"\t\"7\"\t\"8\"\t\"9\"\n"
		"\"0\"\t\"Hello\"\t\"hello\"\t\"hello\"\t\"0\"\t\"0\"\t\"0\"\t\"0\"\t\"0\"\t\"0\"\t\"0\"\t\"0\"\t\"0\"\t\"0\"\n"
		"\"1\"\t\"H\"\t\"e\"\t\"L\"\t\"0.0\"\t\"0.0\"\t\"0.0\"\t\"0.0\"\t\"0.0\"\t\"0.0\"\t\"0.0\"\t\"0.0\"\t\"0.0\"\t\"0.0\"\n"
		"\"2\"\t\"Here's\"\t\"a\"\t\"test\"\t\"100.0\"\t\"100000.254\"\t\"576.22\"\t\"-0.1\"\t\"-100.001\"\t\"100.001\"\t\"1000000\"\t\"4294967296\"\t\"9007199254740992\"\t\"-9007199254740992\"\n"
		"\"3\"\t\"A\"\t\"B\"\t\"C\"\t\"0.00001\"\t\"9999.99999\"\t\"0.00001\"\t\"9999.99999\"\t\"0.00001\"\t\"9999.99999\"\t\"0.00001\"\t\"9999.99999\"\t\"0.00001\"\t\"9999.99999\"\n" ;

		std::size_t number_of_rows = 4 ;
	}
}

namespace {
	void copy_data_to_file( std::string const& data, std::string const& filename ) {
		std::ofstream file( filename.c_str() ) ;
		TEST_ASSERT( file.is_open() ) ;
		file << data ;
		file.close() ;
	}
}

AUTO_TEST_CASE( test_TabDelimitedFormat ) 
{
	std::cerr << "Testing tab-delimited format..." ;
	std::string filename = std::tmpnam(0) ;

	std::cerr << "1..." ;

	{
		copy_data_to_file( globals::data, filename ) ;
		statfile::DelimitedStatSource stat_source_1( filename, "\t" ) ;
		int index ;
		std::string column1, column2, column3 ;
		std::vector< double > doubles( 10 ) ;

		int count = 0 ;
		while( stat_source_1 >> index ) {
			TEST_ASSERT( index == count++ ) ;
			stat_source_1 >> column1 >> column2 >> column3 ;
			for( std::size_t i = 0; i < doubles.size(); ++i ) {
				stat_source_1 >> doubles[i] ;
			}
			stat_source_1 >> statfile::end_row() ;
		}
		
		TEST_ASSERT( stat_source_1.number_of_rows_read() == globals::number_of_rows ) ;
	}
	
	std::cerr << "2..." ;

	{
		copy_data_to_file( globals::quoted_data, filename ) ;
		statfile::DelimitedStatSource stat_source_1( filename, "\t" ) ;
		int index ;
		std::string column1, column2, column3 ;
		std::vector< double > doubles( 10 ) ;

		int count = 0 ;
		while( stat_source_1 >> index ) {
			TEST_ASSERT( index == count++ ) ;
			stat_source_1 >> column1 >> column2 >> column3 ;
			for( std::size_t i = 0; i < doubles.size(); ++i ) {
				stat_source_1 >> doubles[i] ;
			}
			stat_source_1 >> statfile::end_row() ;
		}
		
		TEST_ASSERT( stat_source_1.number_of_rows_read() == globals::number_of_rows ) ;
	}
	
	std::cout << "ok.\n" ;
}

AUTO_TEST_CASE( test_CommaDelimitedFormat ) 
{
	std::cerr << "Testing comma-delimited format..." ;
	
	std::string filename = std::tmpnam(0) ;

//	std::cout << "First filename: \"" << filename << "\".  Second filename: \"" << filename2 << "\".\n" ;

	{
		std::string data = globals::data ;
		std::replace( data.begin(), data.end(), '\t', ',' ) ;
		copy_data_to_file( data, filename ) ;
		statfile::DelimitedStatSource stat_source_1( filename, "," ) ;
		int index ;
		std::string column1, column2, column3 ;
		std::vector< double > doubles( 10 ) ;

		int count = 0 ;
		while( stat_source_1 >> index ) {
			TEST_ASSERT( index == count++ ) ;
			stat_source_1 >> column1 >> column2 >> column3 ;
			for( std::size_t i = 0; i < doubles.size(); ++i ) {
				stat_source_1 >> doubles[i] ;
			}
			stat_source_1 >> statfile::end_row() ;
		}
		
		TEST_ASSERT( stat_source_1.number_of_rows_read() == globals::number_of_rows ) ;
	}

	{
		std::string data = globals::quoted_data ;
		std::replace( data.begin(), data.end(), '\t', ',' ) ;
		copy_data_to_file( data, filename ) ;
		statfile::DelimitedStatSource stat_source_1( filename, "," ) ;
		int index ;
		std::string column1, column2, column3 ;
		std::vector< double > doubles( 10 ) ;

		int count = 0 ;
		while( stat_source_1 >> index ) {
			TEST_ASSERT( index == count++ ) ;
			stat_source_1 >> column1 >> column2 >> column3 ;
			for( std::size_t i = 0; i < doubles.size(); ++i ) {
				stat_source_1 >> doubles[i] ;
			}
			stat_source_1 >> statfile::end_row() ;
		}
		
		TEST_ASSERT( stat_source_1.number_of_rows_read() == globals::number_of_rows ) ;
	}

	std::cout << "ok.\n" ;
}

AUTO_TEST_CASE( test_column_names ) {
		std::cerr << "Testing column names..." ;
		std::string filename = std::tmpnam(0) ;
	//	std::cout << "First filename: \"" << filename << "\".  Second filename: \"" << filename2 << "\".\n" ;
		{
			copy_data_to_file( globals::data, filename ) ;
			statfile::DelimitedStatSource stat_source( filename, "\t" ) ;
			std::vector< std::string > column_names = stat_source.column_names() ;
			TEST_ASSERT( column_names.size() == 14 ) ;
			TEST_ASSERT( column_names[0] == "index" ) ;
			TEST_ASSERT( column_names[1] == "Column1" ) ;
			TEST_ASSERT( column_names[2] == "Column2" ) ;
			TEST_ASSERT( column_names[3] == "Column3" ) ;
			TEST_ASSERT( column_names[4] == "0" ) ;
			TEST_ASSERT( column_names[5] == "1" ) ;
			TEST_ASSERT( column_names[6] == "2" ) ;
			TEST_ASSERT( column_names[7] == "3" ) ;
			TEST_ASSERT( column_names[8] == "4" ) ;
			TEST_ASSERT( column_names[9] == "5" ) ;
			TEST_ASSERT( column_names[10] == "6" ) ;
			TEST_ASSERT( column_names[11] == "7" ) ;
			TEST_ASSERT( column_names[12] == "8" ) ;
			TEST_ASSERT( column_names[13] == "9" ) ;
		}
		
		{
			std::string data = globals::data ;
			std::replace( data.begin(), data.end(), '\t', ',' ) ;
			copy_data_to_file( data, filename ) ;
			statfile::DelimitedStatSource stat_source( filename, "," ) ;
			std::vector< std::string > column_names = stat_source.column_names() ;
			TEST_ASSERT( column_names.size() == 14 ) ;
			TEST_ASSERT( column_names[0] == "index" ) ;
			TEST_ASSERT( column_names[1] == "Column1" ) ;
			TEST_ASSERT( column_names[2] == "Column2" ) ;
			TEST_ASSERT( column_names[3] == "Column3" ) ;
			TEST_ASSERT( column_names[4] == "0" ) ;
			TEST_ASSERT( column_names[5] == "1" ) ;
			TEST_ASSERT( column_names[6] == "2" ) ;
			TEST_ASSERT( column_names[7] == "3" ) ;
			TEST_ASSERT( column_names[8] == "4" ) ;
			TEST_ASSERT( column_names[9] == "5" ) ;
			TEST_ASSERT( column_names[10] == "6" ) ;
			TEST_ASSERT( column_names[11] == "7" ) ;
			TEST_ASSERT( column_names[12] == "8" ) ;
			TEST_ASSERT( column_names[13] == "9" ) ;
		}
		std::cout << "ok.\n" ;
}

