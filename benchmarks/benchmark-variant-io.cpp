#include <boost/variant.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <boost/progress.hpp>

struct MissingValue {} ;

int main() {
	
	std::size_t const N = 100 ;

	std::string number_string ;

	{
		std::stringstream str ;
		std::ifstream f( "/tmp/test/numbers.txt" ) ;

		if( !f.is_open() ) {
			std::cerr << "Please create \"/tmp/test/numbers.txt\".\n" ;
			return -1 ;
		}
		std::size_t count = 0 ;
		std::string entry ;
		while( f >> entry ) {
			str << entry << " " ;
			++count ;
		}
		std::cerr << "Read " << count << " numbers from \"/tmp/test/numbers.txt\".\n" ;
		number_string = str.str() ;
	}
	
	std::cerr << "Numbers begin with: \"" << number_string.substr(0, 50) << "...\".\n" ;

	// Read everything once first to deal with any cache issues.
	{
		std::vector< double > numbers ;
		numbers.reserve( 100000 ) ;
		std::istringstream istr( number_string ) ;
		numbers.clear() ;
		double d ;
		while( istr >> d ) {
			numbers.push_back( d ) ;
		}
	}
	
	std::cerr << "Beginning tests...\n" ;
	
	{
		std::vector< double > numbers ;
		numbers.reserve( 100000 ) ;
		{
			double total_time = 0.0 ;
			for( std::size_t i = 0; i < N; ++i ) {
				numbers.clear() ;
				std::vector< char > v( number_string.begin(), number_string.end() ) ;
				// there's a trailing space in the number string.  Overwrite with a null terminator.
				v.back() = char(0) ;
				boost::timer timer ;
				char* startptr = &v[0] ;
				char* endptr = startptr ;
				for( ; endptr != &(v.back()) ; startptr = endptr ) {
					numbers.push_back( strtod( startptr, &endptr ) ) ;
				}
				total_time += timer.elapsed() ;
			}
			std::cerr << "  strtod reads: " << N << "x" << numbers.size() << " numbers took " << total_time << "s in total, average " << total_time / ( double(N) ) << " numbers per second.\n" ;
		}
	}

	typedef boost::variant< double, int, std::string, MissingValue > Entry ;

	{
		std::vector< Entry > numbers ;
		numbers.reserve( 100000 ) ;
		{
			boost::timer timer ;
			double total_time = 0.0 ;
			for( std::size_t i = 0; i < N; ++i ) {
				std::istringstream istr( number_string ) ;
				numbers.clear() ;
				double d ;
				boost::timer timer ;
				while( istr >> d ) {
					numbers.push_back( d ) ;
				}
				total_time += timer.elapsed() ;
			}
			std::cerr << " Variant reads: " << N << "x" << numbers.size() << " numbers took " << total_time << "s in total, average " << total_time / double(N) << " reads per second.\n" ;
		}
		
		{
			double total_time = 0.0 ;
			for( std::size_t i = 0; i < N; ++i ) {
				std::ostringstream ostr ;
				boost::timer timer ;
				for( std::size_t i = 0; i < numbers.size(); ++i ) {
					ostr << boost::get< double >( numbers[i] ) << " " ;
				}
				total_time += timer.elapsed() ;
			}
			std::cerr << "Variant writes: took " << total_time << "s, average writing time was " << total_time / N << ".\n" ;
		}
		
	}

	{
		std::vector< double > numbers ;
		numbers.reserve( 100000 ) ;
		{
			double total_time = 0.0 ;
			for( std::size_t i = 0; i < N; ++i ) {
				std::istringstream istr( number_string ) ;
				numbers.clear() ;
				double d ;
				boost::timer timer ;
				while( istr >> d ) {
					numbers.push_back( d ) ;
				}
				total_time += timer.elapsed() ;
			}
			std::cerr << "   Plain reads: " << N << "x" << numbers.size() << " numbers took " << total_time << "s in total, average " << total_time / ( double(N) ) << " numbers per second.\n" ;
		}
		
		double total_time = 0.0 ;
		for( std::size_t i = 0; i < N; ++i ) {
			std::ostringstream ostr ;
			boost::timer timer ;
			for( std::size_t i = 0; i < numbers.size(); ++i ) {
				ostr << numbers[i] << " " ;
			}
			total_time += timer.elapsed() ;
		}
		std::cerr << "  Plain writes: took " << total_time << "s, average writing time was " << total_time / N << ".\n" ;
	}

	return 0 ;
}