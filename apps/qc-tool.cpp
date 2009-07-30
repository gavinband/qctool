/*
 * This program, qc-tool, selects rows from a GEN file according to certain criteria.
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
#include <memory>
#include <numeric>
#include "Timer.hpp"
#include "GenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "OptionProcessor.hpp"
#include "RowCondition.hpp"
#include "SNPInListCondition.hpp"
#include "SampleInListCondition.hpp"
#include "Whitespace.hpp"
#include "FileUtil.hpp"
#include "GenRowStatistics.hpp"
#include "SampleRowStatistics.hpp"
#include "GenRowSource.hpp"
#include "ObjectSource.hpp"
#include "GenRowSink.hpp"
#include "SNPDataSource.hpp"
#include "SNPDataSourceChain.hpp"
#include "SNPDataSink.hpp"
#include "SampleInputFile.hpp"
#include "SampleOutputFile.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "wildcard.hpp"
#include "string_utils.hpp"
#include "parse_utils.hpp"
#include "progress_bar.hpp"

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
	    options[ "-g" ]
	        .set_description( "Path of gen file to input" )
	        .set_is_required()
			.set_takes_values()
			.set_maximum_number_of_repeats( 100 )
			.add_value_preprocessor( &expand_filename_wildcards )
	        .add_value_checker( &check_files_are_readable ) ;

	    options[ "-s" ]
	        .set_description( "Path of sample file to input" )
	        .set_takes_single_value()
	        .add_value_checker( &check_files_are_readable ) ;

	    options[ "-og" ]
	        .set_description( "Path of gen file to output" )
	        .set_takes_single_value() ;

		options[ "-os" ]
	        .set_description( "Path of sample file to output" )
	        .set_takes_single_value() ;

	    options[ "-snp-stats" ]
	        .set_description( "Output snp-wise statistics to the given file." )
	        .set_takes_single_value() ;

	    options[ "-sample-stats" ]
	        .set_description( "Output sample-wise statistics to the given file." )
	        .set_takes_single_value() ;

		// SNP filtering options
		options[ "-hwe"]
			.set_description( "Filter out SNPs with HWE exact test statistics less than or equal to the value specified.")
			.set_takes_single_value() ;
		options[ "-snp-missing-rate"]
			.set_description( "Filter out SNPs with missing data rate greater than or equal to the value specified.")
			.set_takes_single_value() ;
		options[ "-snp-interval"]
			.set_description( "Filter out SNPs with position outside the interval [a,b], where a and b are the first and second supplied values" )
			.set_number_of_values_per_use( 2 ) ;
		options[ "-maf"]
			.set_description( "Filter out SNPs whose minor allele frequency lies outside the interval [a,b], where a and b are the first and second supplied values." )
			.set_number_of_values_per_use( 2 ) ;
		options[ "-snp-incl-list"]
			.set_description( "Filter out SNPs whose SNP ID or RSID does not lie in the given file (which must contain a list of whitespace-separated strings)")
			.set_takes_single_value() ;
		options[ "-snp-excl-list"]
			.set_description( "Filter out SNPs whose SNP ID or RSID lies in the given file (which must contain a list of whitespace-separated strings)")
			.set_takes_single_value() ;

		// Sample filtering options
		options[ "-sample-missing-rate" ]
			.set_description( "Filter out samples with missing data rate greater than the value specified.  Note that a full-genome set of GEN files must be supplied.")
			.set_takes_single_value() ;
		options[ "-heterozygosity" ]
			.set_description( "Filter out samples with heterozygosity outside the inteval [a,b], where a and b are the first and second supplied values" )
			.set_number_of_values_per_use( 2 ) ;
		options[ "-sample-incl-list"]
			.set_description( "Filter out samples whose sample ID does not lie in the given file (which must contain a list of whitespace-separated strings)")
			.set_takes_single_value() ;
		options[ "-sample-excl-list"]
			.set_description( "Filter out samples whose sample ID lies in the given file (which must contain a list of whitespace-separated strings)")
			.set_takes_single_value() ;

		options[ "-snp-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in genstat file" )
			.set_takes_single_value()
			.set_default_value( "SNPID, RSID, position, alleles, MAF, HWE, missing" ) ;

		options[ "-sample-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in samplestat file" )
			.set_takes_single_value()
			.set_default_value( std::string("ID1, ID2, missing, heterozygosity") ) ;
			
		options [ "-force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
	}

	GenSelectProcessor( OptionProcessor const& options )
		: m_cout( std::cout.rdbuf() ),
		  m_options( options )
	{
		write_start_banner( m_cout ) ;
		setup() ;
	}

	~GenSelectProcessor()
	{
		write_end_banner( m_cout ) ;
	}
	
private:
	void setup() {
		m_ignore_warnings = m_options.check_if_option_was_supplied( "-force" ) ;
		get_required_filenames() ;
		try {
			open_gen_row_source() ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_cout << "The GEN files specified did not all have the same sample size.\n" ;
			throw ;
		}
		open_sample_row_source() ;
		construct_snp_statistics() ;
		construct_sample_statistics() ;
		construct_snp_filter() ;
		construct_sample_filter() ;
	}

	void get_required_filenames() {
		if( m_options.check_if_option_was_supplied( "-g" ) ) {
			m_gen_filenames = m_options.get_values< std::string >( "-g" ) ;
		}
		if( m_options.check_if_option_was_supplied( "-s" ) ) {
			m_sample_filename = m_options.get_value< std::string >( "-s" ) ;
		}
		if( m_options.check_if_option_was_supplied( "-og" ) ) {
			m_gen_output_filename = m_options.get_value< std::string >( "-og" ) ;
		}
		if( m_options.check_if_option_was_supplied( "-snp-stats" ) ) {
			m_gen_statistic_filename = m_options.get_value< std::string >( "-snp-stats" ) ;
		}
		if( m_options.check_if_option_was_supplied( "-sample-stats" ) ) {
			m_sample_statistic_filename = m_options.get_value< std::string >( "-sample-stats" ) ;
		}
		if( m_options.check_if_option_was_supplied( "-os" )) {
			m_sample_output_filename = m_options.get_value< std::string >( "-os" ) ;
		}
		if( m_options.check_if_option_was_supplied( "-oss" ) ) {
			m_sample_statistic_filename = m_options.get_value< std::string >( "-oss" ) ;
		}
	}

	void open_gen_row_source() {
		Timer timer ;
		
		std::auto_ptr< genfile::SNPDataSourceChain > chain( new genfile::SNPDataSourceChain() ) ;
		m_gen_file_snp_counts.resize( m_gen_filenames.size() ) ;
		for( std::size_t i = 0; i < m_gen_filenames.size(); ++i ) {
		Timer file_timer ;
			m_cout << "(Opening gen file \"" << m_gen_filenames[i] << "\"...)" << std::flush ;
			try {
				std::auto_ptr< genfile::SNPDataSource > snp_data_source( genfile::SNPDataSource::create( m_gen_filenames[i] )) ;
				m_gen_file_snp_counts[i] = snp_data_source->total_number_of_snps() ;
				chain->add_source( snp_data_source ) ;
			}
			catch ( genfile::FileHasTwoConsecutiveNewlinesError const& e ) {
				std::cerr << "\n!!ERROR: a GEN file was specified having two consecutive newlines.\n"
					<< "!! NOTE: popular editors, such as vim and nano, automatically add an extra newline to the file (which you can't see).\n"
					<< "!!     : Please check that each SNP in the file is terminated by a single newline.\n" ;
				throw ;
			}
			m_cout << " (" << file_timer.elapsed() << "s\n" ;
		}

		m_number_of_samples_from_gen_file = chain->number_of_samples() ;
		m_total_number_of_snps = chain->total_number_of_snps() ;

		std::auto_ptr< genfile::SNPDataSource > snp_data_source( chain.release() ) ; 
		gen_row_source.reset( new SNPDataSourceGenRowSource( snp_data_source )) ;

		if( timer.elapsed() > 1.0 ) {
			m_cout << "Opened " << m_gen_filenames.size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
		}
	}

	void open_gen_row_sink() {
		gen_row_sink.reset( new NullObjectSink< GenRow >() ) ;
		if( m_gen_output_filename != "" ) {
			std::auto_ptr< genfile::SNPDataSink > snp_data_sink( genfile::SNPDataSink::create( m_gen_output_filename )) ;
			gen_row_sink.reset( new SNPDataSinkGenRowSink( snp_data_sink )) ;
		}
	}

	void close_gen_row_sink() {
		gen_row_sink.reset( new NullObjectSink< GenRow >() ) ;
	}

	void open_gen_stats_file() {
		// Output GEN stats to std output if no file is supplied.
		// genStatisticOutputFile = OUTPUT_FILE_PTR( new std::ostream( m_cout.rdbuf() )) ;
		if( m_gen_statistic_filename != "" ) {
			genStatisticOutputFile = open_file_for_output( m_gen_statistic_filename ) ;
		}
	}

	void open_sample_row_source() {
		m_sample_row_source.reset( new NullObjectSource< SampleRow >()) ;
		m_have_sample_file = false ;
		if( m_sample_filename != "" ) {
			m_sample_row_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_sample_filename ))) ;
			m_have_sample_file = true ;
		}
	} ;


	void open_sample_row_sink() {
		m_sample_row_sink.reset( new NullObjectSink< SampleRow >() ) ;
		if( m_sample_output_filename != "" ) {
			m_sample_row_sink.reset( new SampleOutputFile< SimpleFileObjectSink< SampleRow > >( open_file_for_output( m_sample_output_filename ))) ;
		}
	}

	void open_sample_stats_file() {
		// Output sample stats to std out if no file is supplied
		// sampleStatisticOutputFile = OUTPUT_FILE_PTR( new std::ostream( std::cout.rdbuf() )) ;
		if( m_sample_statistic_filename != "" ) {
			sampleStatisticOutputFile = open_file_for_output( m_sample_statistic_filename ) ;
		}
	}

	void construct_snp_statistics() {
		std::vector< std::string > row_statistics_specs = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-snp-statistics" ), "," ) ;
		GenRowStatisticFactory::add_statistics( row_statistics_specs, m_row_statistics ) ;
	}

	void construct_sample_statistics() {
		std::vector< std::string > sample_statistics_specs = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-sample-statistics" ), "," ) ;
		SampleRowStatisticFactory::add_statistics( sample_statistics_specs, m_sample_statistics ) ;
	}

	void construct_snp_filter() {
		std::auto_ptr< AndRowCondition > snp_filter( new AndRowCondition() ) ;

		if( m_options.check_if_option_was_supplied( "-hwe" ) ) {
			add_one_arg_condition_to_filter< StatisticGreaterThan >( *snp_filter, "HWE", m_options.get_value< double >( "-hwe" )) ;
		}

		if( m_options.check_if_option_was_supplied( "-snp-missing-rate" ) ) {
			add_one_arg_condition_to_filter< StatisticLessThan >( *snp_filter, "missing", m_options.get_value< double >( "-snp-missing-rate" )) ;
		}

		if( m_options.check_if_option_was_supplied( "-snp-interval" ) ) {
			add_two_arg_condition_to_filter< StatisticInInclusiveRange >( *snp_filter, "snp-position", m_options.get_values< double >( "-snp-interval" )) ;
		}

		if( m_options.check_if_option_was_supplied( "-maf" ) ) {
			add_two_arg_condition_to_filter< StatisticInInclusiveRange >( *snp_filter, "MAF", m_options.get_values< double >( "-maf" )) ;
		}

		if( m_options.check_if_option_was_supplied( "-snp-incl-list" ) ) {
			std::string filename = m_options.get_value< std::string >( "-snp-incl-list" ) ;
			std::auto_ptr< RowCondition > snp_incl_condition( new SNPInListCondition( filename )) ;
			snp_filter->add_subcondition( snp_incl_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "-snp-excl-list" ) ) {
			std::string filename = m_options.get_value< std::string >( "-snp-excl-list" ) ;
			std::auto_ptr< RowCondition > snp_incl_condition( new SNPInListCondition( filename )) ;
			std::auto_ptr< RowCondition > snp_excl_condition( new NotRowCondition( snp_incl_condition )) ;
			snp_filter->add_subcondition( snp_excl_condition ) ;
		}
		
		m_snp_filter = snp_filter ;
	}

	void construct_sample_filter() {
		std::auto_ptr< AndRowCondition > sample_filter( new AndRowCondition() ) ;
		
		if( m_options.check_if_option_was_supplied( "-sample-missing-rate" ) ) {
			add_one_arg_condition_to_filter< StatisticLessThan >( *sample_filter, "missing", m_options.get_value< double >( "-sample-missing-rate" )) ;
		}

		if( m_options.check_if_option_was_supplied( "-heterozygosity" ) ) {
			add_two_arg_condition_to_filter< StatisticInInclusiveRange >( *sample_filter, "heterozygosity", m_options.get_values< double >( "-heterozygosity" )) ;
		}

		if( m_options.check_if_option_was_supplied( "-sample-incl-list" ) ) {
			std::string filename = m_options.get_value< std::string >( "-sample-incl-list" ) ;
			std::auto_ptr< RowCondition > sample_incl_condition( new SampleInListCondition( filename )) ;
			sample_filter->add_subcondition( sample_incl_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "-sample-excl-list" ) ) {
			std::string filename = m_options.get_value< std::string >( "-sample-excl-list" ) ;
			std::auto_ptr< RowCondition > sample_incl_condition( new SampleInListCondition( filename )) ;
			std::auto_ptr< RowCondition > sample_excl_condition( new NotRowCondition( sample_incl_condition )) ;
			sample_filter->add_subcondition( sample_excl_condition ) ;
		}
		
		m_sample_filter = sample_filter ;
	}

	template< typename ConditionType >
	void add_one_arg_condition_to_filter( AndRowCondition& filter, std::string const& statistic_name, double value ) {
		std::auto_ptr< RowCondition > condition( new ConditionType( statistic_name, value )) ;
		filter.add_subcondition( condition ) ;
	}

	template< typename ConditionType >
	void add_two_arg_condition_to_filter( AndRowCondition& filter, std::string const& statistic_name, std::vector< double > values ) {
		assert( values.size() == 2 ) ;
		std::auto_ptr< RowCondition > condition( new ConditionType( statistic_name, values[0], values[1] )) ;
		filter.add_subcondition( condition ) ;
	}

public:
	
	void write_start_banner( std::ostream& oStream ) const {
		oStream << "\nWelcome to qc-tool\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}

	void write_end_banner( std::ostream& oStream ) const {
		oStream << "\nThank you for using qc-tool.\n" ;
	}
	
	void write_preamble( std::ostream& oStream ) const {
		oStream << std::string( 72, '=' ) << "\n\n" ;
		try {
			oStream << std::setw(30) << "Input GEN files:" ;
			for( std::size_t i = 0; i < m_gen_filenames.size(); ++i ) {
				if( i > 0 ) {
					oStream << std::string( 30, ' ' ) ;
				}
				oStream << "  (" << std::setw(6) << m_gen_file_snp_counts[i] << " snps)  " ;
				oStream << "\"" << m_gen_filenames[i] << "\"\n" ;
			}
			oStream << std::string( 30, ' ' ) << "  (total " << m_total_number_of_snps << " snps).\n" ;
			oStream << std::setw(30) << "Output GEN files:"
				<< "  \"" << m_gen_output_filename << "\".\n" ;
			oStream << std::setw(30) << "Input SAMPLE files:"
				<< "  \"" << m_sample_filename << "\".\n" ;
			oStream << std::setw(30) << "Output SAMPLE files:"
				<< "  \"" << m_sample_output_filename << "\".\n" ;
			oStream << std::setw(30) << "SNP statistic output file:"
				<< "  \"" << m_gen_statistic_filename << "\".\n" ;
			oStream << std::setw(30) << "Sample statistic output file:"
				<< "  \"" << m_sample_statistic_filename << "\".\n" ;
			oStream << std::setw(30) << "Sample filter:" 
				<< "  " << *m_sample_filter << ".\n" ;
			oStream << std::setw(30) << "SNP filter:"
				<< "  " << *m_snp_filter << ".\n" ;
			oStream << "\n" ;
		
			if( !m_errors.empty() ) {
				for( std::size_t i = 0; i < m_errors.size(); ++i ) {
					oStream << "!! ERROR: " << m_errors[i] << "\n\n" ;
				}
				oStream << "!! Please correct the above errors and re-run qc-tool.\n\n" ;
				throw GenSelectProcessorException( "Errors were encountered." ) ;
			}


			if( !m_warnings.empty() ) {
				for( std::size_t i = 0; i < m_warnings.size(); ++i ) {
					oStream << "!! WARNING: " << m_warnings[i] << "\n\n" ;
				}
				if( m_ignore_warnings ) {
					oStream << "!! Warnings encountered, but proceeding anyway as -force was supplied.\n\n" ;
				}
				else {
					oStream << "!! To proceed anyway, please run again with the -force option.\n\n" ;
					throw GenSelectProcessorException( "Warnings were encountered." ) ;
				}
			}

			oStream << std::string( 72, '=' ) << "\n\n" ;
		}
		catch (...) {
			oStream << std::string( 72, '=' ) << "\n\n" ;
			throw ;
		}
	}

	void do_checks() {
		for( std::size_t i = 0; i < m_gen_filenames.size(); ++i ) {
			if( strings_are_nonempty_and_equal( m_gen_output_filename, m_gen_filenames[i] )) {
				m_errors.push_back( "Output GEN file \"" + m_gen_output_filename +"\" also specified as input GEN file." ) ;
				break ;
			}
			if( strings_are_nonempty_and_equal( m_sample_output_filename, m_gen_filenames[i] )) {
				m_errors.push_back( "Output SAMPLE file \"" + m_gen_output_filename +"\" also specified as input GEN file." ) ;
				break ;
			}
			if( strings_are_nonempty_and_equal( m_gen_statistic_filename, m_gen_filenames[i] )) {
				m_errors.push_back( "Output GEN statistic file \"" + m_gen_output_filename +"\" also specified as input GEN file." ) ;
				break ;
			}
			if( strings_are_nonempty_and_equal( m_sample_statistic_filename, m_gen_filenames[i] )) {
				m_errors.push_back( "Output SAMPLE statistic file \"" + m_gen_output_filename +"\" also specified as input GEN file." ) ;
				break ;
			}
		}
		if( strings_are_nonempty_and_equal( m_sample_output_filename, m_sample_filename )) {
			m_errors.push_back( "Output sample file \"" + m_sample_filename + "\" also specified as input sample file." ) ;
		}
		if( strings_are_nonempty_and_equal( m_sample_output_filename, m_gen_output_filename )) {
			m_errors.push_back( "The GEN and SAMPLE output filenames must differ." ) ;
		}
		if( strings_are_nonempty_and_equal( m_gen_statistic_filename, m_sample_statistic_filename )) {
			m_errors.push_back( "The gen statistic and sample statistic filenames must differ." ) ;
		}
		if( m_sample_output_filename != "" && m_gen_filenames.size() != 23 ) {
			m_warnings.push_back( "You are outputting a sample file, but the number of gen files is not 23.") ;
		}
		if( m_sample_statistic_filename != "" && m_gen_filenames.size() != 23 ) {
			m_warnings.push_back( "You are outputting a sample statistic file, but the number of gen files is not 23.") ;
		}
		if( m_sample_statistic_filename != "" && m_sample_filename == "" ) {
			m_warnings.push_back( "You are outputting a sample statistic file, but no input sample file has been supplied.") ;
		}
		if( m_gen_output_filename == "" && m_sample_output_filename == "" && m_gen_statistic_filename == "" && m_sample_statistic_filename == "" ) {
			m_warnings.push_back( "You have not specified any output files.  This will produce only console output." ) ;
		}
		if( (m_snp_filter->number_of_subconditions() == 0) && (m_sample_filter->number_of_subconditions() == 0) && (m_gen_statistic_filename == "") && (m_sample_statistic_filename == "") ) {
			m_warnings.push_back(
			"You have not specified any filters, nor any statistic output files.\n"
			"This will just copy input to output files (taking into account any\n"
			"file format changes).  You can do this faster with gen-convert." ) ;
		}
	}
	
	bool strings_are_nonempty_and_equal( std::string const& left, std::string const& right ) {
		return (!left.empty()) && (!right.empty()) && (left == right) ;
	}
	
	void process() {
		try {
			unsafe_process() ;
		}
		catch( StatisticNotFoundException const& e ) {
			std::cerr << "!! ERROR: " << e << ".\n" ;
			std::cerr << "Note: required statistics must be added using -statistics.\n" ;
		}
		catch( GToolException const& e) {
			std::cerr << "!! ERROR: " << e << ".\n" ;
		}
	}

private:

	void unsafe_process() {
		filter_sample_rows() ;
		process_gen_rows() ;
		process_sample_rows() ;
	}

	void filter_sample_rows() {
		m_cout << "Processing sample file...\n" ;
		open_sample_row_source() ;
		SampleRow sample_row ;
		m_number_of_sample_file_rows = 0 ;
		for( ; (*m_sample_row_source) >> sample_row; ++m_number_of_sample_file_rows ) {
			if( m_sample_filter->check_if_satisfied( sample_row )) {
				m_indices_of_filtered_in_sample_rows.push_back( m_number_of_sample_file_rows ) ;
			}
			else {
				m_indices_of_filtered_out_sample_rows.push_back( m_number_of_sample_file_rows ) ;
			}
		}
		
		m_cout << "Total samples " << m_number_of_sample_file_rows << ", keeping " << m_indices_of_filtered_in_sample_rows.size() << ", filtering out " << m_indices_of_filtered_out_sample_rows.size() << ".\n" ;
	}

	void process_gen_rows() {
		m_cout << "Processing GEN file(s)...\n" ;
		Timer timer ;
		open_gen_row_sink() ;
		open_gen_stats_file() ;

		if( genStatisticOutputFile.get() ) { 
			*genStatisticOutputFile << "        " ;
			m_row_statistics.format_column_headers( *genStatisticOutputFile ) << "\n";
		}

		double last_time = -1 ;
		InternalStorageGenRow row ;
		row.set_number_of_samples( m_number_of_samples_from_gen_file ) ;
		std::size_t number_of_snps_processed = 0 ;
		while( read_gen_row( row )) {
			preprocess_gen_row( row ) ;
			process_gen_row( row, ++number_of_snps_processed ) ;
			accumulate_per_column_amounts( row, m_per_column_amounts ) ;
			double time_now = timer.elapsed() ;
			if( (time_now - last_time > 1.0) || (number_of_snps_processed == m_total_number_of_snps) ) {
				m_cout
					<< "\r"
					<< get_progress_bar( 30, static_cast< double >( number_of_snps_processed ) / m_total_number_of_snps )
					<< " (" << number_of_snps_processed << " / " << m_total_number_of_snps
					<< ", " << std::fixed << std::setprecision(1) << time_now << "s)"
					<< std::flush ;
				last_time = time_now ;
			}
		}
		m_cout << "\n" ;
		assert( number_of_snps_processed = m_total_number_of_snps ) ;
		std::cerr << "\nProcessed GEN file(s) (" << m_total_number_of_snps << " rows) in " << timer.elapsed() << " seconds.\n" ;
	
		m_cout << "Post-processing (updating file header, compression)..." << std::flush ;
		timer.restart() ;
		// Close the output gen file(s) now
		close_gen_row_sink() ;
		std::cerr << " (" << timer.elapsed() << "s)\n" ;
	}

	bool read_gen_row( GenRow& row ) {
		return (*gen_row_source) >> row ;
	}

	void preprocess_gen_row( InternalStorageGenRow& row ) const {
		check_gen_row( row ) ;
		row.filter_out_samples_with_indices( m_indices_of_filtered_out_sample_rows ) ;
	}
	
	void check_gen_row( GenRow& row ) const {
		check_gen_row_has_correct_number_of_samples( row ) ;
	}

	void check_gen_row_has_correct_number_of_samples( GenRow& row ) const {
		if( m_have_sample_file ) {
			if( row.number_of_samples() != m_number_of_sample_file_rows ) {
				throw GenAndSampleFileMismatchException( "GEN file and sample file have mismatching number of samples." ) ;
			}
		}
	}

	void process_sample_rows() {
		Timer timer ;
		
		// re-open sample row source.
		open_sample_row_source() ;
		open_sample_row_sink() ;
		open_sample_stats_file() ;
		
		if( m_have_sample_file ) {
			if( m_number_of_sample_file_rows != m_per_column_amounts.size() ) {
				throw GenAndSampleFileMismatchException( "Sample file and GEN file have mismatching number of samples." ) ;
			}
		}
		
		if( sampleStatisticOutputFile.get() ) {
			*sampleStatisticOutputFile << "        " ;
			m_sample_statistics.format_column_headers( *sampleStatisticOutputFile ) << "\n" ;
		}

		SampleRow sample_row ;

		std::size_t i = 0 ;
		for( ; i < m_per_column_amounts.size(); ++i ) {
			if( m_have_sample_file ) {
				(*m_sample_row_source) >> sample_row ;
				assert( *m_sample_row_source ) ; // assume sample file has not changed since count_sample_rows()
			}

			m_sample_statistics.process( sample_row, m_per_column_amounts[i], m_total_number_of_snps ) ;

			if( sampleStatisticOutputFile.get() ) {
				*sampleStatisticOutputFile << std::setw(8) << (i+1) << m_sample_statistics << "\n" ;
			}
			m_sample_statistics.add_to_sample_row( sample_row, "missing" ) ;
			m_sample_statistics.add_to_sample_row( sample_row, "heterozygosity" ) ;
			(*m_sample_row_sink) << sample_row ;
		}

		std::cerr << "Processed sample file (" << i << " rows) in " << timer.elapsed() << " seconds.\n" ;
	}

	void process_gen_row( GenRow const& row, std::size_t row_number ) {
		m_row_statistics.process( row ) ;
		m_row_statistics.round_genotype_amounts() ;
		if( m_snp_filter->check_if_satisfied( m_row_statistics )) {
			output_gen_row( row ) ;
			output_gen_row_stats( m_row_statistics, row_number ) ;
		}
	}

	void output_gen_row( GenRow const& row ) {
		(*gen_row_sink) << row ;
	}

	void output_gen_row_stats( GenotypeAssayStatistics const& row_statistics, std::size_t row_number ) {
		if( genStatisticOutputFile.get() ) {
			if( m_row_statistics.size() > 0 ) {
				(*genStatisticOutputFile)
					<< std::setw(8) << std::left << row_number
					<< m_row_statistics << "\n" ;
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
	
private:
	
	std::ostream m_cout ;
	
	std::auto_ptr< ObjectSource< GenRow > > gen_row_source ;
	std::auto_ptr< ObjectSink< GenRow > > gen_row_sink ;
	std::auto_ptr< ObjectSource< SampleRow > > m_sample_row_source ;
	std::auto_ptr< ObjectSink< SampleRow > > m_sample_row_sink ;
	OUTPUT_FILE_PTR genStatisticOutputFile ;
	OUTPUT_FILE_PTR sampleStatisticOutputFile ;
	OptionProcessor const& m_options ;
	
	GenRowStatistics m_row_statistics ;
	SampleRowStatistics m_sample_statistics ;
	std::auto_ptr< AndRowCondition > m_snp_filter ;
	std::auto_ptr< AndRowCondition > m_sample_filter ;
	
	std::vector< std::size_t > m_gen_file_snp_counts ;
	std::size_t m_total_number_of_snps ;
	std::size_t m_number_of_samples_from_gen_file ;
	std::size_t m_number_of_sample_file_rows ;
	bool m_have_sample_file ;
	
	std::vector< GenotypeProportions > m_per_column_amounts ;

	std::vector< std::size_t > m_indices_of_filtered_in_sample_rows ;
	std::vector< std::size_t > m_indices_of_filtered_out_sample_rows ;
	
	std::vector< std::string > m_gen_filenames ;
	std::string m_sample_filename ;
	std::string m_gen_output_filename ;
	std::string m_sample_output_filename ;
	std::string m_gen_statistic_filename ;
	std::string m_sample_statistic_filename ;
	
	std::vector< std::string > m_warnings ;
	std::vector< std::string > m_errors ;
	bool m_ignore_warnings ;
} ;

int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		GenSelectProcessor::declare_options( options ) ;
		options.process( argc, argv ) ;
    }
    catch( std::exception const& exception ) {
        std::cerr << "!! Error: " << exception.what() << ".\n";
        std::cerr << "Usage: qc-tool [options]\n"
                << options
                << "\n" ;
        return -1 ;
    }

	// main section

	try	{
		GenSelectProcessor processor( options ) ;
		processor.do_checks() ;
		processor.write_preamble( std::cout ) ;
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
	return parse_integer_in_half_open_range( a_string, 1, 101 ) != 101 ;
}


std::vector< std::string > expand_filename_wildcards( std::string const& option_name, std::vector< std::string > const& filenames ) {
	std::vector< std::string > result ;
	for( std::size_t i = 0; i < filenames.size(); ++i ) {
		if( filenames[i].find( '#' ) != std::string::npos ) {
			std::pair< std::vector< std::string >, std::vector< std::string > > expanded_filename = find_files_matching_path_with_wildcard( filenames[i], '#' ) ;
			if( expanded_filename.first.empty() ) {
				throw OptionValueInvalidException( option_name, filenames, "No file can be found matching filename \"" + filenames[i] + "\"." ) ;
			}

			// We only allow matches corresponding to numbers 1 to 100
			for( std::size_t j = 0; j < expanded_filename.first.size(); ) {
				if( !check_if_string_is_a_number_from_1_to_100( expanded_filename.second[j] )) {
					expanded_filename.first.erase( expanded_filename.first.begin() + j ) ;
					expanded_filename.second.erase( expanded_filename.second.begin() + j ) ;
				}
				else {
					++j ;
				}
			}

			std::copy( expanded_filename.first.begin(), expanded_filename.first.end(), std::back_inserter( result )) ;
		}
		else {
			result.push_back( filenames[i] ) ;
		}
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



