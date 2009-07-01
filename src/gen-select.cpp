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
#include "SNPInListCondition.hpp"
#include "Whitespace.hpp"
#include "FileUtil.hpp"
#include "GenRowStatistics.hpp"
#include "SampleRowStatistics.hpp"
#include "GenRowFileSource.hpp"
#include "GenRowFileSink.hpp"
#include "SampleInputFile.hpp"
#include "SampleOutputFile.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "HardyWeinbergExactTestStatistic.hpp"
#include "string_utils.hpp"

std::vector< std::string > expand_filename_wildcards( std::string const& option_name, std::vector< std::string > const& filenames ) ;
void check_files_are_readable( std::string const& option_name, std::vector< std::string > const& filenames ) ;
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
		
		// File options		
	    options[ "--g" ]
	        .set_description( "Path of gen file to input" )
	        .set_is_required()
			.set_takes_values()
			.set_maximum_number_of_repeats( 100 )
			.add_value_preprocessor( &expand_filename_wildcards )
	        .add_value_checker( &check_files_are_readable ) ;

	    options[ "--s" ]
	        .set_description( "Path of sample file to input" )
	        .set_takes_single_value()
	        .add_value_checker( &check_files_are_readable ) ;

	    options[ "--og" ]
	        .set_description( "Path of gen file to output" )
	        .set_takes_single_value() ;

	    options[ "--ogs" ]
	        .set_description( "Path of gen statistic file to output" )
	        .set_takes_single_value() ;

		options[ "--os" ]
	        .set_description( "Path of sample file to output" )
	        .set_takes_single_value() ;

		// Statistics-related options

		options[ "--snp-stats" ]
	        .set_description( "Output per-snp statistics." ) ;
		options[ "--sample-stats" ]
	        .set_description( "Output per-sample statistics for missing rate and heterozygosity" ) ;

		// SNP filtering options
		options[ "--hwe"]
			.set_description( "Filter out SNPs with HWE exact test statistics less than or equal to the value specified.")
			.set_takes_single_value() ;
		options[ "--snp-missing-rate"]
			.set_description( "Filter out SNPs with missing data rate greater than or equal to the value specified.")
			.set_takes_single_value() ;
		options[ "--snp-interval"]
			.set_description( "Filter out SNPs with position outside the interval [a,b], where a and b are the first and second supplied values" )
			.set_number_of_values_per_use( 2 ) ;
		options[ "--maf"]
			.set_description( "Filter out SNPs whose minor allele frequency lies outside the interval [a,b], where a and b are the first and second supplied values." )
			.set_number_of_values_per_use( 2 ) ;
		options[ "--snp-incl-list"]
			.set_description( "Filter out SNPs whose SNP ID or RSID does not lie in the given file (which must contain a list of whitespace-separated strings)") ;
		options[ "--snp-excl-list"]
			.set_description( "Filter out SNPs whose SNP ID or RSID lies in the given file (which must contain a list of whitespace-separated strings)") ;

		// Sample filtering options
		options[ "--sample-missing-rate" ]
			.set_description( "Filter out samples with missing data rate greater than the value specified.  Note that a full-genome set of GEN files must be supplied.") ;
		options[ "--heterozygosity" ]
			.set_description( "Filter out samples with heterozygosity outside the inteval [a,b], where a and b are the first and second supplied values" )
			.set_number_of_values_per_use( 2 ) ;
		options[ "--sample-incl-list"]
			.set_description( "Filter out samples whose sample ID does not lie in the given file (which must contain a list of whitespace-separated strings)") ;
		options[ "--sample-excl-list"]
			.set_description( "Filter out samples whose sample ID lies in the given file (which must contain a list of whitespace-separated strings)") ;

		options[ "--row-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in genstat file" )
			.set_takes_single_value()
			.set_default_value( "SNPID, RSID, position, alleles, MAF, HWE, missing-rate" ) ;

		options[ "--sample-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in samplestat file" )
			.set_takes_single_value()
			.set_default_value( std::string("missing-rate, heterozygosity") ) ;
	}

	GenSelectProcessor( OptionProcessor const& options )
		: m_options( options )
	{
		setup() ;
	}
	
	void open_sample_row_source() {
		m_sample_row_source.reset( new NullObjectSource< SampleRow >()) ;
		m_have_sample_file = false ;
		if( m_options.check_if_option_was_supplied( "--s" )) {
			std::string sampleFileName = m_options.get_value< std::string >( "--s" ) ;
			m_sample_row_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( sampleFileName ))) ;
			m_have_sample_file = true ;
		}
	} ;

	void open_gen_row_source() {
		std::vector< std::string > gen_filenames = m_options.get_values< std::string >( "--g" ) ;
		gen_row_source = get_genrow_source_from_files( gen_filenames ) ;
		std::cout << "gen-select: opened the following GEN files for input:\n" ;
		for( std::size_t i = 0; i < gen_filenames.size(); ++i ) {
			std::cout << "  - \"" << gen_filenames[i] << "\".\n" ;
		}
	}
	
	void open_gen_row_sink() {
		gen_row_sink.reset( new NullObjectSink< GenRow >() ) ;
		if( m_options.check_if_option_was_supplied( "--og" ) ) {
			std::string genOutputFileName = m_options.get_value< std::string >( "--og" ) ;
			if( determine_file_mode( genOutputFileName ) & e_BinaryMode ) {
				gen_row_sink.reset( new SimpleGenRowBinaryFileSink( open_file_for_output( genOutputFileName ))) ;
			} else {
				gen_row_sink.reset( new SimpleFileObjectSink< GenRow >( open_file_for_output( genOutputFileName ))) ;
			}
		}
	}

	void open_sample_row_sink() {
		m_sample_row_sink.reset( new NullObjectSink< SampleRow >() ) ;
		if( m_options.check_if_option_was_supplied( "--os" )) {
			std::string sampleOutputFileName = m_options.get_value< std::string >( "--os" ) ;
			m_sample_row_sink.reset( new SampleOutputFile< SimpleFileObjectSink< SampleRow > >( open_file_for_output( sampleOutputFileName ))) ;
		}
		else {
			OUTPUT_FILE_PTR sampleOutputFile( new std::ostream( std::cout.rdbuf() ) ) ;
			m_sample_row_sink.reset( new SampleOutputFile< SimpleFileObjectSink< SampleRow > >( sampleOutputFile )) ;
		}
	}

	void construct_snp_filter() {
		std::map< std::string, std::string > filter_components ;
		std::auto_ptr< AndRowCondition > snp_filter( new AndRowCondition() ) ;

		if( m_options.check_if_option_was_supplied( "--hwe" ) ) {
			double threshhold = m_options.get_value< double >( "--hwe" ) ;
			std::auto_ptr< RowCondition > hwe_condition( new GenotypeAssayStatisticGreaterThan( "HWE", threshhold )) ;
			snp_filter->add_subcondition( hwe_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "--snp-missing-rate" ) ) {
			double threshhold = m_options.get_value< double >( "--snp-missing-rate" ) ;
			std::auto_ptr< RowCondition > missing_rate_condition( new GenotypeAssayStatisticLessThan( "missing-rate", threshhold )) ;
			snp_filter->add_subcondition( missing_rate_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "--snp-interval" ) ) {
			std::vector< double > range = m_options.get_values< double >( "--snp-interval" ) ;
			assert( range.size() == 2 ) ;
			assert( range[0] <= range[1] ) ;
			std::auto_ptr< RowCondition > position_condition( new GenotypeAssayStatisticInInclusiveRange( "snp-position", range[0], range[1] )) ;
			snp_filter->add_subcondition( position_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "--maf" ) ) {
			std::vector< double > range = m_options.get_values< double >( "--maf" ) ;
			assert( range.size() == 2 ) ;
			assert( range[0] <= range[1] ) ;
			std::auto_ptr< RowCondition > maf_condition( new GenotypeAssayStatisticInInclusiveRange( "MAF", range[0], range[1] )) ;
			snp_filter->add_subcondition( maf_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "--snp-incl-list" ) ) {
			std::string filename = m_options.get_value< std::string >( "--snp-incl-list" ) ;
			std::auto_ptr< RowCondition > snp_incl_condition( new SNPInListCondition( filename )) ;
			snp_filter->add_subcondition( snp_incl_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "--snp-excl-list" ) ) {
			std::string filename = m_options.get_value< std::string >( "--snp-excl-list" ) ;
			std::auto_ptr< RowCondition > snp_incl_condition( new SNPInListCondition( filename )) ;
			std::auto_ptr< RowCondition > snp_excl_condition( new NotRowCondition( snp_incl_condition )) ;
			snp_filter->add_subcondition( snp_excl_condition ) ;
		}
		
		m_snp_filter = snp_filter ;
		
		std::cout << "gen-select: I will keep SNPs which satisfy: "
			<< (*m_snp_filter) << ".\n" ;
	}

	void construct_sample_filter() {
		if( m_options.check_if_option_was_supplied( "--snp-incl-list" ) ) {
			std::string filename = m_options.get_value< std::string >( "--snp-incl-list" ) ;
			std::auto_ptr< SampleCondition > snp_incl_condition( new SampleInListCondition( filename )) ;
			snp_filter->add_subcondition( snp_incl_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "--snp-excl-list" ) ) {
			std::string filename = m_options.get_value< std::string >( "--snp-excl-list" ) ;
			std::auto_ptr< SampleCondition > snp_incl_condition( new SampleInListCondition( filename )) ;
			std::auto_ptr< SampleCondition > snp_excl_condition( new NotSampleCondition( snp_incl_condition )) ;
			snp_filter->add_subcondition( snp_excl_condition ) ;
		}
	}

	void setup() {
		open_gen_row_source() ;
		open_gen_row_sink() ;
		open_sample_row_source() ;
		open_sample_row_sink() ;
		
		// Output Gen row stats to cout by default.
		genStatisticOutputFile = OUTPUT_FILE_PTR( new std::ostream( std::cout.rdbuf() )) ;
		sampleStatisticOutputFile = OUTPUT_FILE_PTR( new std::ostream( std::cout.rdbuf() )) ;

		if( m_options.check_if_option_was_supplied( "--ogs" ) ) {
			std::string genStatisticFileName = m_options.get_value< std::string >( "--ogs" ) ;
			genStatisticOutputFile = open_file_for_output( genStatisticFileName ) ;
		}
		
		std::vector< std::string > row_statistics_specs = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "--row-statistics" ), "," ) ;
		GenRowStatisticFactory::add_statistics( row_statistics_specs, m_row_statistics ) ;

		std::vector< std::string > sample_statistics_specs = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "--sample-statistics" ), "," ) ;
		SampleRowStatisticFactory::add_statistics( sample_statistics_specs, m_sample_statistics ) ;
		
		// Get the condition to apply to rows.
		construct_snp_filter() ;
	}
	
	void process() {
		try {
			unsafe_process() ;
		}
		catch( StatisticNotFoundException const& e ) {
			std::cerr << "!! ERROR: " << e << ".\n" ;
			std::cerr << "Note: required statistics must be added using --statistics.\n" ;
		}
	}

private:

	void unsafe_process() {
		count_sample_rows() ;
		process_gen_rows() ;
		process_sample_rows() ;
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
		
		while( (*gen_row_source) >> row ) {
			check_gen_row( row ) ;
			process_gen_row( row, ++m_number_of_gen_rows ) ;
			if( m_number_of_gen_rows % 1000 == 0 ) {
				std::cout << "Processed " << m_number_of_gen_rows << " rows (" << timer.elapsed() << "s)...\n" ;
			}
			accumulate_per_column_amounts( row, m_per_column_amounts ) ;
		}
	#ifdef HAVE_BOOST_TIMER
		std::cerr << "gen-select: processed GEN file(s) (" << m_number_of_gen_rows << " rows) in " << timer.elapsed() << " seconds.\n" ;
	#endif
	}

	void check_gen_row( GenRow& row ) {
		// no-op
		if( m_have_sample_file ) {
			if( row.number_of_samples() != m_number_of_sample_file_rows ) {
				throw GenAndSampleFileMismatchException( "GEN file and sample file have mismatching number of samples." ) ;
			}
		}
	}

	void process_gen_row( GenRow& row, std::size_t row_number ) {
		m_row_statistics.process( row ) ;
		m_row_statistics.round_genotype_amounts() ;
		if( m_snp_filter->check_if_satisfied( row, &m_row_statistics )) {
			(*gen_row_sink) << row ;

			if( genStatisticOutputFile.get() ) {
				if( m_row_statistics.size() > 0 ) {
					(*genStatisticOutputFile)
						<< std::setw(8) << std::left << row_number
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
	
	void count_sample_rows() {
		open_sample_row_source() ;
		SampleRow sample_row ;
		m_number_of_sample_file_rows = 0 ;
		while( (*m_sample_row_source) >> sample_row ) {
			++m_number_of_sample_file_rows ;
		}
	}

	void process_sample_rows() {
		#ifdef HAVE_BOOST_TIMER
			boost::timer timer ;
		#endif
		
		// re-open sample row source.
		open_sample_row_source() ;
		
		if( m_have_sample_file ) {
			if( m_number_of_sample_file_rows != m_per_column_amounts.size() ) {
				throw GenAndSampleFileMismatchException( "Sample file and GEN file have mismatching number of samples." ) ;
			}
		}
		
		*sampleStatisticOutputFile << "        " ;
		m_sample_statistics.format_column_headers( *sampleStatisticOutputFile ) << "\n" ;

		SampleRow sample_row ;

		std::size_t i = 0 ;
		for( ; i < m_per_column_amounts.size(); ++i ) {
			if( m_have_sample_file ) {
				(*m_sample_row_source) >> sample_row ;
				assert( *m_sample_row_source ) ; // assume sample file has not changed since count_sample_rows()
			}

			m_sample_statistics.process( sample_row, m_per_column_amounts[i], m_number_of_gen_rows ) ;
			*sampleStatisticOutputFile << std::setw(8) << (i+1) << m_sample_statistics << "\n" ;
			m_sample_statistics.add_to_sample_row( sample_row ) ;
			(*m_sample_row_sink) << sample_row ;
		}

		#ifdef HAVE_BOOST_TIMER
			std::cerr << "gen-select: processed sample file (" << i << " rows) in " << timer.elapsed() << " seconds.\n" ;
		#endif
	}

private:
	
	std::auto_ptr< ObjectSource< GenRow > > gen_row_source ;
	std::auto_ptr< ObjectSink< GenRow > > gen_row_sink ;
	std::auto_ptr< ObjectSource< SampleRow > > m_sample_row_source ;
	std::auto_ptr< ObjectSink< SampleRow > > m_sample_row_sink ;
	OUTPUT_FILE_PTR genStatisticOutputFile ;
	OUTPUT_FILE_PTR sampleStatisticOutputFile ;
	OptionProcessor const& m_options ;
	
	GenRowStatistics m_row_statistics ;
	SampleRowStatistics m_sample_statistics ;
	std::auto_ptr< RowCondition > m_snp_filter ;
	
	std::size_t m_number_of_gen_rows ;
	std::size_t m_number_of_sample_file_rows ;
	bool m_have_sample_file ;
	
	std::vector< GenotypeProportions > m_per_column_amounts ;
} ;

int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		GenSelectProcessor::declare_options( options ) ;
		options.process( argc, argv ) ;
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
		GenSelectProcessor processor( options ) ;
		processor.process() ;
	}
	catch( std::exception const& e )
	{
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}


bool check_if_string_is_a_number_from_1_to_100( std::string const& a_string ) {
	std::istringstream inStream( a_string ) ;
	int i ;
	inStream >> i ;
	inStream.peek() ;
	if( !inStream.eof()) {
		return false ;
	}
	if( i < 1 || i > 100 ) {
		return false ;
	}
	return true ;
}


std::vector< std::string > expand_filename_wildcards( std::string const& option_name, std::vector< std::string > const& filenames ) {
	std::vector< std::string > result ;
	for( std::size_t i = 0; i < filenames.size(); ++i ) {
		std::vector< std::string > expanded_filename = find_files_matching_path_with_wildcard( filenames[i], '#', &check_if_string_is_a_number_from_1_to_100 ) ;
		if( expanded_filename.empty() ) {
			throw OptionValueInvalidException( option_name, filenames, "No file can be found matching filename \"" + filenames[i] + "\"." ) ;
		}
		std::copy( expanded_filename.begin(), expanded_filename.end(), std::back_inserter( result )) ;
	}

	return result ;
}

void check_files_are_readable( std::string const& option_name, std::vector< std::string > const& filenames ) {
	if( filenames.empty() ) {
		throw OptionValueInvalidException( option_name, filenames, "At least one filename must be supplied.\"" ) ;
	}
	for( std::size_t i = 0; i < filenames.size(); ++i ) {
    	std::ifstream file( filenames[i].c_str() ) ;
	    if( !file.good() ) {
	        throw OptionValueInvalidException( option_name, filenames, "File \"" + filenames[i] + "\" is not readable." ) ;
	    }
	}
}

void check_condition_spec( std::string const& option_name, std::vector< std::string > const& condition_specs ) {
	for( std::size_t i = 0; i < condition_specs.size() ; ++i ) {
		if( condition_specs[i].size() == 0 ) {
			throw OptionValueInvalidException( option_name, condition_specs, "Condition spec \"" + condition_specs[i] + "\" supplied for option " + option_name + " must be nonempty." ) ;
		}
	}
}



