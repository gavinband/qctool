/*
 * This program, sample-select, selects columns from a GEN file (and samples from the corresponding sample file)
 * according to certain criteria.
 *
 * Program arguments:
 *
 */

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <numeric>
#include "GenRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "OptionProcessor.hpp"
#include "RowCondition.hpp"
#include "RowConditionFactory.hpp"
#include "Whitespace.hpp"
#include "FileUtil.hpp"
#include "GenRowStatistics.hpp"
#include "GenRowSimpleFileSource.hpp"
#include "GenRowSimpleFileSink.hpp"
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
        .set_takes_single_value()
        .set_value_checker( &check_if_file_is_readable ) ;

    options[ "--s" ]
        .set_description( "Path of sample file to input" )
        .set_takes_single_value()
        .set_value_checker( &check_if_file_is_readable ) ;

    options[ "--og" ]
        .set_description( "Path of gen file to output" )
        .set_takes_single_value() ;

    options[ "--ogs" ]
        .set_description( "Path of gen statistic file to output" )
        .set_takes_single_value() ;

		options[ "--so" ]
        .set_description( "Path of sample file to output" )
        .set_takes_single_value() ;

	options[ "--samples" ]
        .set_description( "Sample selector" ) ;

	options[ "--condition" ]
        .set_description( "Condition spec for gen file row selection" )
		.set_takes_single_value()
		.set_default_value( std::string("") ) ;

	options[ "--statistics" ]
        .set_description( "Comma-seperated list of statistics to calculate in genstat file" )
		.set_takes_single_value()
		.set_default_value( std::string("") ) ;

		options.process( argc, argv ) ;
}        


void process_gen_rows( INPUT_FILE_PTR genFile, OUTPUT_FILE_PTR genOutputFile, OUTPUT_FILE_PTR genStatisticOutputFile, OptionProcessor const& options ) {
	Whitespace whitespace ;
	GenRow row ;
	GenRowStatistics row_statistics ;
	
	{
		std::string statistics_spec_str = options.get_argument_value< std::string >( "--statistics" ) ;
		std::vector< std::string > statistics_specs = split( statistics_spec_str, "," ) ;
		for( std::size_t i = 0; i < statistics_specs.size(); ++i ) {
			std::string spec = strip( statistics_specs[i] ) ;
			if( spec != "" ) {
				row_statistics.add_statistic( spec, GenotypeAssayStatisticFactory::create_statistic( spec )) ;
			}
		}
	}

	if( genStatisticOutputFile.get() ) { 
		(*genStatisticOutputFile) << "SNPID          RSID        A B : " ;
		row_statistics.format_column_headers( *genStatisticOutputFile ) << "\n";
	}

	// Get the condition to apply to rows.
	std::string condition_spec = options.get_argument_value< std::string >( "--condition" ) ;
	std::auto_ptr< RowCondition > condition = RowConditionFactory::create_condition( condition_spec ) ;

	try {		
		GenRowSimpleFileSource row_source( genFile ) ;
		GenRowSimpleFileSink row_sink( genFile ) ;

		while( !row_source.check_if_empty() ) {
			row_source >> row ;
			row_statistics.process( row ) ;
			row_statistics.round_genotype_amounts() ;

			if( condition->check_if_satisfied( row, &row_statistics )) {
				row_sink << row ;

				if( genStatisticOutputFile.get() ) {
					(*genStatisticOutputFile)
						<< std::setw(14) << std::left << row.SNPID() << " "
						<< std::setw(11) << std::left << row.RSID() << " "
						<< row.first_allele() << " " << row.second_allele() << " : "
						<< row_statistics << "\n" ;
				}
			}
		}

		if( genStatisticOutputFile.get() ) { 
			(*genStatisticOutputFile) << "SNPID          RSID        A B : " ;
			row_statistics.format_column_headers( *genStatisticOutputFile ) << "\n";
		}
	}
	catch( StatisticNotFoundException const& e ) {
		std::cerr << "!! ERROR: " << e << ".\n" ;
		std::cerr << "Note: statistics required by the --condition argument must be added using --statistics.\n" ;

	}
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
			std::string genFileName = options.get_argument_value< std::string >( "--g" ) ;
			INPUT_FILE_PTR genFile = open_file_for_input( genFileName.c_str() ) ;
			OUTPUT_FILE_PTR genOutputFile, genStatisticOutputFile ;
			if( options.check_if_option_was_supplied( "--og" ) ) {
				std::string genOutputFileName = options.get_argument_value< std::string >( "--og" ) ;
				genOutputFile = open_file_for_output( genOutputFileName.c_str() ) ;
			}
			if( options.check_if_option_was_supplied( "--ogs" ) ) {
				std::string genStatisticFileName = options.get_argument_value< std::string >( "--ogs" ) ;
				genStatisticOutputFile = open_file_for_output( genStatisticFileName.c_str() ) ;
			}
			else {
				genStatisticOutputFile = OUTPUT_FILE_PTR( new std::ostream( std::cout.rdbuf() )) ;
			}

			process_gen_rows( genFile, genOutputFile, genStatisticOutputFile, options ) ;
		}
	}
	catch( GToolException const& e )
	{
		std::cerr << "!! Error: " << e.message() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}
 

