/*
 * This program, bechmark-statistics, reads some GEN rows from a file an
 *
 * Program arguments:
 *
 */

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <numeric>
#include "../config.hpp"
#if HAVE_BOOST_TIMER
	#include <boost/timer.hpp>
	typedef boost::timer Timer ;
#else
	struct Timer
	{
		double elapsed() const { return 0.0 ;}
	} ;
#endif
#include "GenRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "OptionProcessor.hpp"
#include "FileUtil.hpp"
#include "GenRowStatistics.hpp"
#include "SNPDataSource.hpp"
#include "GenRowSource.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "HardyWeinbergExactTestStatistic.hpp"
#include "string_utils.hpp"

void process_options( OptionProcessor& options, int argc, char** argv ) {
    options[ "--g" ]
        .set_description( "Path of gen file to input" )
        .set_is_required()
		.set_takes_single_value() ;

	options[ "--nr" ]
        .set_description( "Number of rows to read" )
		.set_takes_single_value()
		.set_default_value( 10000 ) ;

	options.process( argc, argv ) ;
}        


void print_results( std::ostream&, std::map< double, std::vector< std::string > > const& results ) {
	// Write output
	std::cout << "--- RESULTS ---\n" ;
	std::map< double, std::vector< std::string > >::const_iterator
		i( results.begin() ),
		end_i( results.end() ) ;

	for( ; i != end_i; ++i ) {
		std::vector< std::string > const& stats = i->second ;
		for( std::size_t j = 0; j < stats.size(); ++j )
			std::cout << std::setw(25) << std::left << stats[j] + ":" << std::fixed << std::setprecision(5) << i->first << "\n" ;
	}
}


void process_gen_rows( ObjectSource< GenRow >& gen_row_source, OptionProcessor const& options ) {
	GenRow row ;

	std::size_t number_of_rows_to_read = options.get_value< std::size_t >( "--nr" ) ;

	std::cout << "Reading " << number_of_rows_to_read << " GenRows...\n" ;
	std::size_t count = 0;
	Timer timer ;
	// Read in a list of genRows.
	while(( count < number_of_rows_to_read ) && ( gen_row_source >> row )) {
		++count ;
	}
	
	std::cout << "Read " << count << " rows in " << timer.elapsed() << "s.\n" ;
}


int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
        process_options( options, argc, argv ) ;
    }
    catch( std::exception const& exception ) {
        std::cerr << "!! Error: " << exception.what() << ".\n";
        std::cerr << "Usage: gen-select [options]\n"
                << options
                << "\n" ;
        return -1 ;
    }

	// main section

	try	{
		Whitespace whitespace ;

		{
			std::auto_ptr< ObjectSource< GenRow > > genRowSource( new NullObjectSource< GenRow >() ) ;

			if( options.check_if_option_was_supplied( "--g" )) {
				std::string genFileName = options.get_value< std::string >( "--g" ) ;
				genRowSource.reset( new SNPDataSourceGenRowSource( genfile::SNPDataSource::create( genFileName ))) ;
			}

			process_gen_rows( *genRowSource, options ) ;
		}
	}
	catch( GToolException const& e )
	{
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}
 

