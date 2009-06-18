/*
 * This program, gen-select, selects rows from a GEN file according to certain criteria.
 * - rows where some genotype data is missing, or none is.
 * - rows where genotype data is, or is not in hardy-weinberg equilibrium.
 * - rows with SNP IDs in a given list.
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
#endif
#include "GenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "OptionProcessor.hpp"
#include "RowCondition.hpp"
#include "RowConditionFactory.hpp"
#include "Whitespace.hpp"
#include "FileUtil.hpp"
#include "GenRowStatistics.hpp"
#include "SimpleFileObjectSource.hpp"
#include "SimpleFileObjectSink.hpp"
#include "SampleInputFile.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "HardyWeinbergExactTestStatistic.hpp"
#include "string_utils.hpp"

void check_if_file_is_readable( std::string const& option_name, std::string const& filename ) {
    std::ifstream file( filename.c_str() ) ;
    if( !file.good() ) {
        throw ArgumentInvalidException( "File \"" + filename + "\" supplied for option " + option_name + " is not readable." ) ;
    }    
}

void check_condition_spec( std::string const& option_name, std::string const& condition_spec ) {
	if( condition_spec.size() == 0 ) {
		throw ArgumentInvalidException( "Condition spec \"" + condition_spec + "\" supplied for option " + option_name + " must be nonempty." ) ;
	}
}

void process_options( OptionProcessor& options, int argc, char** argv ) {
    options[ "--g" ]
        .set_description( "Path of gen file to input" )
        .set_is_required()
        .set_takes_value()
        .set_value_checker( &check_if_file_is_readable ) ;

    options[ "--s" ]
        .set_description( "Path of sample file to input" )
        .set_takes_value()
        .set_value_checker( &check_if_file_is_readable ) ;

    options[ "--og" ]
        .set_description( "Path of gen file to output" )
        .set_takes_value() ;

    options[ "--ogs" ]
        .set_description( "Path of gen statistic file to output" )
        .set_takes_value() ;

	options[ "--so" ]
        .set_description( "Path of sample file to output" )
        .set_takes_value() ;

	options[ "--samples" ]
        .set_description( "Sample selector" ) ;

	options[ "--condition" ]
        .set_description( "Condition spec for gen file row selection" )
		.set_takes_value()
		.set_default_value( std::string("") ) ;

	options[ "--statistics" ]
        .set_description( "Comma-seperated list of statistics to calculate in genstat file" )
		.set_takes_value()
		.set_default_value( std::string("") ) ;

	options.process( argc, argv ) ;
}        


void process_sample_rows( ObjectSource< SampleRow >& sample_row_source, ObjectSink< SampleRow >& sample_row_sink, OptionProcessor const& options ) {
	SampleRow sample_row ;
	while( !sample_row_source.check_if_empty() ) {
		sample_row_source >> sample_row ;
		std::cout << sample_row ;
	}
}

void process_gen_rows( ObjectSource< GenRow >& gen_row_source, ObjectSink< GenRow >& gen_row_sink, OUTPUT_FILE_PTR genStatisticOutputFile, OptionProcessor const& options ) {
#ifdef HAVE_BOOST_TIMER
	boost::timer timer ;
#endif
	Whitespace whitespace ;
	GenRow row ;
	GenRowStatistics row_statistics ;
	
	{
		std::string statistics_spec_str = options.get_value< std::string >( "--statistics" ) ;
		std::vector< std::string > statistics_specs = split( statistics_spec_str, "," ) ;
		for( std::size_t i = 0; i < statistics_specs.size(); ++i ) {
			std::string spec = strip( statistics_specs[i] ) ;
			if( spec != "" ) {
				row_statistics.add_statistic( spec, GenotypeAssayStatisticFactory::create_statistic( spec )) ;
			}
		}
	}

	if( genStatisticOutputFile.get() ) { 
		row_statistics.format_column_headers( *genStatisticOutputFile ) << "\n";
	}

	// Get the condition to apply to rows.
	std::string condition_spec = options.get_value< std::string >( "--condition" ) ;
	std::auto_ptr< RowCondition > condition = RowConditionFactory::create_condition( condition_spec ) ; 	

	try {		
		while( !gen_row_source.check_if_empty() ) {
			gen_row_source >> row ;
			row_statistics.process( row ) ;
			row_statistics.round_genotype_amounts() ;

			if( condition->check_if_satisfied( row, &row_statistics )) {
				gen_row_sink << row ;

				if( genStatisticOutputFile.get() ) {
					(*genStatisticOutputFile)
						<< row_statistics << "\n" ;
				}
			}
		}

		if( genStatisticOutputFile.get() ) { 
			row_statistics.format_column_headers( *genStatisticOutputFile ) << "\n";
		}
	}
	catch( StatisticNotFoundException const& e ) {
		std::cerr << "!! ERROR: " << e << ".\n" ;
		std::cerr << "Note: statistics required by the --condition argument must be added using --statistics.\n" ;

	}
#ifdef HAVE_BOOST_TIMER
	std::cerr << "gen-select: process GEN file in " << timer.elapsed() << " seconds.\n" ;
#endif
}



int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
        process_options( options, argc, argv ) ;
    }
    catch( ArgumentProcessingException const& exception ) {
        std::cerr << "!! Error: " << exception.message() << ".\n";
        std::cerr << "Usage: gen-select [options]\n"
                << options
                << "\n" ;
        return -1 ;
    }

	// main section

	try	{
		Whitespace whitespace ;
		int number_of_snps = 0 ;

		{
			OUTPUT_FILE_PTR genStatisticOutputFile ;
			std::auto_ptr< ObjectSource< GenRow > > genRowSource( new NullObjectSource< GenRow >() ) ;
			std::auto_ptr< ObjectSink< GenRow > > genRowSink( new NullObjectSink< GenRow >() ) ;
			std::auto_ptr< ObjectSource< SampleRow > > sampleRowSource( new NullObjectSource< SampleRow >() ) ;
			std::auto_ptr< ObjectSink< SampleRow > > sampleRowSink( new NullObjectSink< SampleRow >() ) ;

			if( options.check_if_argument_was_supplied( "--g" )) {
				std::string genFileName = options.get_value< std::string >( "--g" ) ;
				genRowSource.reset( new SimpleFileObjectSource< GenRow >( open_file_for_input( genFileName ))) ;
			}
			if( options.check_if_argument_was_supplied( "--s" )) {
				std::string sampleFileName = options.get_value< std::string >( "--s" ) ;
				sampleRowSource.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( sampleFileName ))) ;
			}
			if( options.check_if_argument_was_supplied( "--og" ) ) {
				std::string genOutputFileName = options.get_value< std::string >( "--og" ) ;
				genRowSink.reset( new SimpleFileObjectSink< GenRow >( open_file_for_output( genOutputFileName ))) ;
			}
			if( options.check_if_argument_was_supplied( "--os" ) ) {
				std::string sampleOutputFileName = options.get_value< std::string >( "--os" ) ;
				sampleRowSink.reset( new SimpleFileObjectSink< SampleRow >( open_file_for_output( sampleOutputFileName ))) ;
			}
			if( options.check_if_argument_was_supplied( "--ogs" ) ) {
				std::string genStatisticFileName = options.get_value< std::string >( "--ogs" ) ;
				genStatisticOutputFile = open_file_for_output( genStatisticFileName ) ;
			}
			else {
				// output to std::cout.
				genStatisticOutputFile = OUTPUT_FILE_PTR( new std::ostream( std::cout.rdbuf() )) ;
			}

			process_sample_rows( *sampleRowSource, *sampleRowSink, options ) ;
			process_gen_rows( *genRowSource, *genRowSink, genStatisticOutputFile, options ) ;

		}
	}
	catch( GToolException const& e )
	{
		std::cerr << "!! Error: " << e.message() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}
 

