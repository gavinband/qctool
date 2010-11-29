#include <iostream>
#include <utility>
#include <map>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "test_case.hpp"
#include "statfile/TabDelimitedStatSource.hpp"

namespace globals {
	std::string data =
	"# This is a file\n"
	"# These lines are comments, which should be ignored.\n"
	"#\n"
	"index\tColumn1\tColumn2\tColumn3\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\n"
	"0\tHello\thello\thello\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
	"1\tH\te\tL\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n"
	"100\tHere's\ta\ttest\t100.0\t100000.254\t576.22\t-0.1\t-100.001\t100.001\t1000000\t4294967296\t9007199254740992\t-9007199254740992\n"
	"110\tA\tB\tC\t0.00001\t9999.99999\t0.00001\t9999.99999\t0.00001\t9999.99999\t0.00001\t9999.99999\t0.00001\t9999.99999\n" ;

	int number_of_data_columns = 10 ;
	std::size_t number_of_rows = 4 ;
}

void copy_data_to_file( std::string const& data, std::string const& filename ) {
	std::ofstream file( filename.c_str() ) ;
	TEST_ASSERT( file.is_open() ) ;
	file << data ;
	file.close() ;
}

AUTO_TEST_CASE( test_TabDelimitedFormat ) 
{
	std::string filename = std::tmpnam(0) ;

//	std::cout << "First filename: \"" << filename << "\".  Second filename: \"" << filename2 << "\".\n" ;

	copy_data_to_file( globals::data, filename ) ;

	{
		statfile::TabDelimitedStatSource stat_source_1( filename ) ;
		int index ;
		std::string column1, column2, column3 ;
		std::vector< double > doubles( 10 ) ;

		while( stat_source_1 >> index ) {
			stat_source_1 >> column1 >> column2 >> column3 ;
			for( std::size_t i = 0; i < doubles.size(); ++i ) {
				stat_source_1 >> doubles[i] ;
			}
			stat_source_1 >> statfile::end_row() ;
		}
		
		TEST_ASSERT( stat_source_1.number_of_rows_read() == globals::number_of_rows ) ;
	}

	std::cout << "Success.\n" ;
}

AUTO_TEST_MAIN
{
	test_TabDelimitedFormat() ;
}
