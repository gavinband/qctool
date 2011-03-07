#include <boost/variant.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <boost/progress.hpp>

struct MissingValue {} ;

int main() {
	
	std::size_t const N = 100 ;

	std::stringstream str ;

	{
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
		
	}
	
	std::cerr << "Numbers begin with: \"" << str.str().substr(0, 50) << "...\".\n" ;

	// Read everything once first to deal with any cache issues.
	{
		std::vector< double > numbers ;
		numbers.reserve( 100000 ) ;
		std::istringstream istr( str.str() ) ;
		numbers.clear() ;
		double d ;
		while( istr >> d ) {
			numbers.push_back( d ) ;
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
				std::istringstream istr( str.str() ) ;
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
				std::istringstream istr( str.str() ) ;
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