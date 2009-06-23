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
#include "SampleRowStatistics.hpp"
#include "GenRowFileSource.hpp"
#include "GenRowFileSink.hpp"
#include "SampleInputFile.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "HardyWeinbergExactTestStatistic.hpp"
#include "string_utils.hpp"

void check_if_file_is_readable( std::string const& option_name, std::string const& filename ) ;
void check_condition_spec( std::string const& option_name, std::string const& condition_spec ) ;


struct GenSelectProcessorException: public GToolException
{
	GenSelectProcessorException( std::string const& msg )
		: GToolException( msg )
	{}
} ;

struct GenAndSampleFileMismatchException: public GenSelectProcessorException
{
	GenAndSampleFileMismatchException( std::string const& msg )
		: GenSelectProcessorException( msg )
	{}
} ;


struct GenSelectProcessor
{
public:
	static void declare_options( OptionProcessor & options ) {
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

		options[ "--row-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in genstat file" )
			.set_takes_value()
			.set_default_value( std::string("") ) ;

		options[ "--sample-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in samplestat file" )
			.set_takes_value()
			.set_default_value( std::string("") ) ;
	}

	GenSelectProcessor( OptionProcessor const& options )
		: m_options( options )
	{
		setup() ;
	}
	
	
	
	void open_sample_row_source() {
		m_sample_row_source.reset( new NullObjectSource< SampleRow >()) ;
		if( m_options.check_if_argument_was_supplied( "--s" )) {
			std::string sampleFileName = m_options.get_value< std::string >( "--s" ) ;
			m_sample_row_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( sampleFileName ))) ;
		}
	} ;

	void open_gen_row_source() {
		gen_row_source.reset( new NullObjectSource< GenRow >() ) ;
		std::string genFileName = m_options.get_value< std::string >( "--g" ) ;
		if( determine_file_mode( genFileName ) & e_BinaryMode ) {
			gen_row_source.reset( new SimpleGenRowBinaryFileSource( open_file_for_input( genFileName ))) ;
		} else {
			gen_row_source.reset( new SimpleFileObjectSource< GenRow >( open_file_for_input( genFileName ))) ;
		}
	}
	
	void open_gen_row_sink() {
		gen_row_sink.reset( new NullObjectSink< GenRow >() ) ;
		if( m_options.check_if_argument_was_supplied( "--og" ) ) {
			std::string genOutputFileName = m_options.get_value< std::string >( "--og" ) ;
			if( determine_file_mode( genOutputFileName ) & e_BinaryMode ) {
				gen_row_sink.reset( new SimpleGenRowBinaryFileSink( open_file_for_output( genOutputFileName ))) ;
			} else {
				gen_row_sink.reset( new SimpleFileObjectSink< GenRow >( open_file_for_output( genOutputFileName ))) ;
			}
		}
	}

	void open_sample_row_sink() {
		sample_row_sink.reset( new NullObjectSink< SampleRow >() ) ;
	}

	void setup() {
		open_gen_row_source() ;
		open_gen_row_sink() ;
		open_sample_row_source() ;
		open_sample_row_sink() ;
		
		// Output Gen row stats to cout by default.
		genStatisticOutputFile = OUTPUT_FILE_PTR( new std::ostream( std::cout.rdbuf() )) ;
		sampleStatisticOutputFile = OUTPUT_FILE_PTR( new std::ostream( std::cout.rdbuf() )) ;

		if( m_options.check_if_argument_was_supplied( "--ogs" ) ) {
			std::string genStatisticFileName = m_options.get_value< std::string >( "--ogs" ) ;
			genStatisticOutputFile = open_file_for_output( genStatisticFileName ) ;
		}
		
		std::vector< std::string > row_statistics_specs = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "--row-statistics" ), "," ) ;
		GenRowStatisticFactory::add_statistics( row_statistics_specs, m_row_statistics ) ;

		std::vector< std::string > sample_statistics_specs = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "--sample-statistics" ), "," ) ;
		SampleRowStatisticFactory::add_statistics( sample_statistics_specs, m_sample_statistics ) ;
		
		// Get the condition to apply to rows.
		std::string condition_spec = m_options.get_value< std::string >( "--condition" ) ;
		m_row_condition = RowConditionFactory::create_condition( condition_spec ) ; 	
	}
	
	void process() {
		try {
			process_unsafe() ;
		}
		catch( StatisticNotFoundException const& e ) {
			std::cerr << "!! ERROR: " << e << ".\n" ;
			std::cerr << "Note: statistics required by the --condition argument must be added using --statistics.\n" ;
		}
	}

private:

	void process_unsafe() {
		process_sample_rows() ;
		process_gen_rows() ;
		process_per_sample_statistics() ;
	}

	void process_gen_rows() {
	#ifdef HAVE_BOOST_TIMER
		boost::timer timer ;
	#endif

		if( genStatisticOutputFile.get() ) { 
			*sampleStatisticOutputFile << "        " ;
			m_row_statistics.format_column_headers( *genStatisticOutputFile ) << "\n";
		}

		GenRow row ;
		m_number_of_gen_rows = 0 ;
		
		while( !gen_row_source->check_if_empty() ) {
			(*gen_row_source) >> row ;
			process_gen_row( row, ++m_number_of_gen_rows ) ;
			accumulate_per_column_amounts( row, m_per_column_amounts ) ;
		}
	#ifdef HAVE_BOOST_TIMER
		std::cerr << "gen-select: processed GEN file (" << m_number_of_gen_rows << " rows) in " << timer.elapsed() << " seconds.\n" ;
	#endif
	}
	
	void process_gen_row( GenRow& row, std::size_t row_number ) {
		m_row_statistics.process( row ) ;
		m_row_statistics.round_genotype_amounts() ;
		if( m_row_condition->check_if_satisfied( row, &m_row_statistics )) {
			(*gen_row_sink) << row ;

			if( genStatisticOutputFile.get() ) {
				if( m_row_statistics.size() > 0 ) {
					(*genStatisticOutputFile)
						<< std::setw(8) << std::left << ++row_number
						<< m_row_statistics << "\n" ;
				}
			}
		}
	}

	void accumulate_per_column_amounts( GenRow& row, std::vector< GenotypeProportions >& per_column_amounts ) {
		// Keep totals for per-column stats.
		if( per_column_amounts.empty() ) {
			per_column_amounts.reserve( row.number_of_samples() ) ;
			std::copy( row.begin_genotype_proportions(), row.end_genotype_proportions(), std::back_inserter( per_column_amounts )) ;
		}
		else {
			assert( per_column_amounts.size() == row.number_of_samples() ) ;
			std::transform( per_column_amounts.begin(), per_column_amounts.end(),
			 				row.begin_genotype_proportions(),
			 				per_column_amounts.begin(),
							std::plus< GenotypeProportions >() ) ;
		}
	}
	
	void process_sample_rows() {
		SampleRow sample_row ;
		while( !m_sample_row_source->check_if_empty() ) {
			(*m_sample_row_source) >> sample_row ;
			std::cout << sample_row ;
		}
	}

	void process_per_sample_statistics() {
		#ifdef HAVE_BOOST_TIMER
			boost::timer timer ;
		#endif
		
		// re-open sample row source.
		open_sample_row_source() ;
		SampleRow sample_row ;
		bool mismatch = false ;

		*sampleStatisticOutputFile << "        " ;
		m_sample_statistics.format_column_headers( *sampleStatisticOutputFile ) << "\n" ;
		std::size_t i = 0 ;
		for( ; i < m_per_column_amounts.size(); ++i ) {
			if( !m_sample_row_source->check_if_empty() ) {
				(*m_sample_row_source) >> sample_row ;
			}
			else {
				// We'll throw an exception, but first print out all the stats as this might be useful.
				mismatch = true ;
			}
			m_sample_statistics.process( sample_row, m_per_column_amounts[i], m_number_of_gen_rows ) ;
	
			if( m_sample_statistics.size() > 0 ) {
				(*sampleStatisticOutputFile) << std::setw(8) << std::left << (i+1)
					<< m_sample_statistics << "\n";
			}
		}
		if( mismatch ) {
			throw GenAndSampleFileMismatchException( "The sample file had fewer rows than the gen file has columns (or it didn't exist)." ) ;
		}

		#ifdef HAVE_BOOST_TIMER
			std::cerr << "gen-select: processed sample file (" << i << " rows) in " << timer.elapsed() << " seconds.\n" ;
		#endif
	}

private:
	
	std::auto_ptr< ObjectSource< GenRow > > gen_row_source ;
	std::auto_ptr< ObjectSink< GenRow > > gen_row_sink ;
	std::auto_ptr< ObjectSource< SampleRow > > m_sample_row_source ;
	std::auto_ptr< ObjectSink< SampleRow > > sample_row_sink ;
	OUTPUT_FILE_PTR genStatisticOutputFile ;
	OUTPUT_FILE_PTR sampleStatisticOutputFile ;
	OptionProcessor const& m_options ;
	
	GenRowStatistics m_row_statistics ;
	SampleRowStatistics m_sample_statistics ;
	std::auto_ptr< RowCondition > m_row_condition ;
	
	std::size_t m_number_of_gen_rows ;
	
	std::vector< GenotypeProportions > m_per_column_amounts ;
} ;

int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		GenSelectProcessor::declare_options( options ) ;
		options.process( argc, argv ) ;
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
		GenSelectProcessor processor( options ) ;
		processor.process() ;
	}
	catch( GToolException const& e )
	{
		std::cerr << "!! Error: " << e.message() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}

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



