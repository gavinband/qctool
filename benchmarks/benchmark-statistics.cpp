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
#include "appcontext/OptionProcessor.hpp"
#include "FileUtil.hpp"
#include "GenRowStatistics.hpp"
#include "SimpleFileObjectSource.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "HardyWeinbergExactTestStatistic.hpp"
#include "string_utils/string_utils.hpp"
#include "appcontext/appcontext.hpp"

using namespace appcontext ;

void process_options( appcontext::OptionProcessor& options, int argc, char** argv ) {
    options[ "--g" ]
        .set_description( "Path of gen file to input" )
        .set_is_required()
		.set_takes_single_value() ;

	options[ "--nr" ]
        .set_description( "Number of rows to read" )
		.set_default_value( 100 ) ;

	options[ "--ni" ]
        .set_description( "Number of iterations per statistic" )
		.set_default_value( 100 ) ;

	options[ "--statistics" ]
        .set_description( "Comma-seperated list of statistics to calculate in genstat file" )
		.set_takes_single_value()
		.set_default_value( std::string("") ) ;

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
	// Begin by reading 100 rows from the input file (or as many as possible)
	Whitespace whitespace ;
	InternalStorageGenRow row ;

	std::size_t number_of_rows_to_read = options.get_value< std::size_t >( "--nr" ) ;
	std::size_t number_of_iterations = options.get_value< std::size_t >( "--ni" ) ;

	std::string statistics_spec_str = options.get_value< std::string >( "--statistics" ) ;
	std::vector< std::string > statistic_specs = string_utils::split_and_strip_discarding_empty_entries( statistics_spec_str, "," ) ;

	try {		
		// Read in a list of genRows.
		std::vector<InternalStorageGenRow> genrows ;
		InternalStorageGenRow genrow ;
		{
			std::size_t count = 0;
			while(( count < number_of_rows_to_read ) && ( gen_row_source >> static_cast< GenRow& >(genrow) )) {
				genrows.push_back( genrow ) ;
				++count ;
			}
		}

		// Order results by time, increasing.
		std::map< double, std::vector< std::string > > results ;

		// Do a benchmark for each stat.
		for( std::size_t i = 0; i < statistic_specs.size(); ++i ) {
			std::cout << "Benchmarking statistic \"" << statistic_specs[i] << "\".\n" ;
			GenRowStatistics row_statistics ;
			row_statistics.add_statistic( statistic_specs[i], GenotypeAssayStatisticFactory::create_statistic( statistic_specs[i] )) ;
			appcontext::Timer timer ;
			for( std::size_t n = 0; n < number_of_iterations; ++n ) {
				for( std::size_t j = 0; j < genrows.size(); ++j ) {
					row_statistics.process( genrows[j] ) ;
					row_statistics.round_genotype_amounts() ;
					row_statistics.get_value< double >( statistic_specs[i] ) ;
				}
			}
			results[ timer.elapsed() ].push_back( statistic_specs[i] ) ;
		}

		print_results( std::cout, results) ;
	}
	catch( StatisticNotFoundException const& e ) {
		std::cerr << "!! ERROR: " << e << ".\n" ;
		std::cerr << "Note: statistics required by the --condition argument must be added using --statistics.\n" ;

	}
}


int main( int argc, char** argv ) {
	appcontext::OptionProcessor options ;
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
			OUTPUT_FILE_PTR genStatisticOutputFile ;
			std::auto_ptr< ObjectSource< GenRow > > genRowSource( new NullObjectSource< GenRow >() ) ;

			if( options.check_if_option_was_supplied( "--g" )) {
				std::string genFileName = options.get_value< std::string >( "--g" ) ;
				genRowSource.reset( new SimpleFileObjectSource< GenRow >( open_file_for_input( genFileName ))) ;
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
 

