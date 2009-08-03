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
#include "VerboseMessageWriter.hpp"

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

struct NumberOfSamplesMismatchException: public GenSelectProcessorException
{
	NumberOfSamplesMismatchException()
		: GenSelectProcessorException( "NumberOfSamplesMismatchException" )
	{}
} ;

struct SampleFileHasTooFewSamplesException: public GenSelectProcessorException
{
	SampleFileHasTooFewSamplesException()
		: GenSelectProcessorException( "SampleFileHasTooFewSamplesException" )
	{}
} ;

struct ProblemsWereEncountered: public GenSelectProcessorException
{
	ProblemsWereEncountered()
		: GenSelectProcessorException( "ProblemsWereEncountered" )
	{}
} ;


struct GenSelectProcessor
{
public:
	static void declare_options( OptionProcessor & options ) {

		// Meta-options
		options.set_help_option( "-help" ) ;
		
		// File options	
		options.declare_group( "Data file options" ) ;
	    options[ "-g" ]
	        .set_description( "Path of gen file to input.  Repeat this option, or use the numerical wildcard character '#', to specify several files." )
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

		// Statistic file options
		options.declare_group( "Statistic file options" ) ;
	    options[ "-snp-stats" ]
	        .set_description( "Output snp-wise statistics to the given file." )
	        .set_takes_single_value() ;

	    options[ "-sample-stats" ]
	        .set_description( "Output sample-wise statistics to the given file." )
	        .set_takes_single_value() ;

		options[ "-snp-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in genstat file." )
			.set_takes_single_value()
			.set_default_value( "SNPID, RSID, position, minor_allele, major_allele, MAF, HWE, missing" ) ;

		options[ "-sample-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in samplestat file." )
			.set_takes_single_value()
			.set_default_value( std::string("ID1, ID2, missing, heterozygosity") ) ;

		// SNP filtering options
		options.declare_group( "SNP filtering options" ) ;
		options[ "-hwe"]
			.set_description( "Filter out SNPs with HWE p-value less than or equal to the value specified.")
			.set_takes_single_value() ;
		options[ "-snp-missing-rate"]
			.set_description( "Filter out SNPs with missing data rate greater than or equal to the value specified.")
			.set_takes_single_value() ;
		options[ "-snp-interval"]
			.set_description( "Filter out SNPs with position outside the interval [a,b]." )
			.set_number_of_values_per_use( 2 ) ;
		options[ "-maf"]
			.set_description( "Filter out SNPs whose minor allele frequency lies outside the interval [a,b]." )
			.set_number_of_values_per_use( 2 ) ;
		options[ "-snp-incl-list"]
			.set_description( "Filter out SNPs whose SNP ID or RSID does not lie in the given file.")
			.set_takes_single_value() ;
		options[ "-snp-excl-list"]
			.set_description( "Filter out SNPs whose SNP ID or RSID lies in the given file.")
			.set_takes_single_value() ;

		// Sample filtering options
		options.declare_group( "Sample filtering options" ) ;
		options[ "-sample-missing-rate" ]
			.set_description( "Filter out samples with missing data rate greater than the value specified.")
			.set_takes_single_value() ;
		options[ "-heterozygosity" ]
			.set_description( "Filter out samples with heterozygosity outside the inteval [a,b]." )
			.set_number_of_values_per_use( 2 ) ;
		options[ "-sample-incl-list"]
			.set_description( "Filter out samples whose sample ID does not lie in the given file.")
			.set_takes_single_value() ;
		options[ "-sample-excl-list"]
			.set_description( "Filter out samples whose sample ID lies in the given file.")
			.set_takes_single_value() ;

		options.declare_group( "Other options" ) ;
		options [ "-force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
		options [ "-diagnose-sample-filter" ]
			.set_description( "Print diagnostic information about each filtered out sample." ) ;
		options [ "-diagnose-snp-filter" ]
			.set_description( "Print diagnostic information about each filtered out snp." ) ;

		options.option_excludes_group( "-snp-stats", "SNP filtering options" ) ;
		options.option_excludes_group( "-sample-stats", "Sample filtering options" ) ;
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
		process_other_options() ;
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
		else if( m_options.check_if_option_was_supplied( "-sample-stats" )) {
			// sample output file defaults to overwriting sample file.
			m_sample_output_filename = m_sample_filename ;
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
			m_cout << " (" << file_timer.elapsed() << "s)\n" ;
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
		if( m_sample_filename != "" ) {
			m_sample_row_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_sample_filename ))) ;
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

	void process_other_options() {
		m_ignore_warnings = m_options.check_if_option_was_supplied( "-force" ) ;
		m_diagnose_sample_filter = m_options.check_if_option_was_supplied( "-diagnose-sample-filter" ) ;
		m_diagnose_snp_filter = m_options.check_if_option_was_supplied( "-diagnose-snp-filter" ) ;
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

			oStream << std::setw(30) << "# of samples in input files:"
				<< "  " << m_number_of_samples_from_gen_file << ".\n" ;
			oStream << std::setw(30) << "# of samples after filtering:"
				<< "  " << m_number_of_samples_from_gen_file - m_indices_of_filtered_out_samples.size()
				<< " (" << m_indices_of_filtered_out_samples.size()
				<< " filtered out).\n" ;
			
			oStream << "\n" << std::string( 72, '=' ) << "\n\n" ;

			if( !m_errors.empty() ) {
				for( std::size_t i = 0; i < m_errors.size(); ++i ) {
					oStream << "!! ERROR: " << m_errors[i] << "\n\n" ;
				}
				oStream << "!! Please correct the above errors and re-run qc-tool.\n" ;
				throw ProblemsWereEncountered() ;
			}


			if( !m_warnings.empty() ) {
				for( std::size_t i = 0; i < m_warnings.size(); ++i ) {
					oStream << "!! WARNING: " << m_warnings[i] << "\n\n" ;
				}
				if( m_ignore_warnings ) {
					oStream << "!! Warnings were encountered, but proceeding anyway as -force was supplied.\n" ;
					oStream << "\n" << std::string( 72, '=' ) << "\n\n" ;
				}
				else {
					oStream << "!! Warnings were encountered.  To proceed anyway, please run again with the -force option.\n" ;
					throw ProblemsWereEncountered() ;
				}
			}
		}
		catch (...) {
			throw ;
		}
	}

	void preprocess() {
		preprocess_sample_rows() ;
		check_for_errors_and_warnings() ;
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

	void preprocess_sample_rows() {
		try {
			unsafe_preprocess_sample_rows() ;
		}
		catch( ConditionValueNotFoundException const& ) {
			m_cout << "\n\n!! ERROR: The input sample file must contain entries for all values used to filter on.\n"
				<< "!! This includes \"missing\" and \"heterozygosity\".\n" ;
			throw ;
		}
	}

	void unsafe_preprocess_sample_rows() {
		if( m_sample_filename != "" ) {
			Timer timer ;
			m_cout << "(Preprocessing sample file...)" ;
			open_sample_row_source() ;
			SampleRow sample_row ;
			m_number_of_sample_file_rows = 0 ;
			for( ; (*m_sample_row_source) >> sample_row; ++m_number_of_sample_file_rows ) {
				m_sample_rows.push_back( sample_row ) ;
				if( !m_sample_filter->check_if_satisfied( sample_row )) {
					m_indices_of_filtered_out_samples.push_back( m_number_of_sample_file_rows ) ;
				}
			}

			m_cout << "(" << std::fixed << std::setprecision(1) << timer.elapsed() << "s)\n" ;
			if( m_number_of_sample_file_rows != m_number_of_samples_from_gen_file ) {
				throw NumberOfSamplesMismatchException() ;
			}
		}
		else {
			m_sample_rows.resize( m_number_of_samples_from_gen_file ) ;
		}
	}

	void check_for_errors_and_warnings() {
		check_for_errors() ;
		check_for_warnings() ;
	}

	void check_for_errors() {
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
		if( strings_are_nonempty_and_equal( m_sample_output_filename, m_gen_output_filename )) {
			m_errors.push_back( "The GEN and SAMPLE output filenames must differ." ) ;
		}
		if( strings_are_nonempty_and_equal( m_gen_statistic_filename, m_sample_statistic_filename )) {
			m_errors.push_back( "The gen statistic and sample statistic filenames must differ." ) ;
		}
		if( m_sample_filename == "" && m_sample_filter->number_of_subconditions() != 0 ) {
			m_errors.push_back( "To filter on samples, please supply an input sample file." ) ;
		}
	}
	
	void check_for_warnings() {
		if( (m_sample_output_filename != "" || m_sample_statistic_filename != "") && m_gen_filenames.size() != 23 ) {
			m_warnings.push_back( "You are outputting a sample or sample statistic file, but the number of gen files is not 23.\n"
			"   (I suspect there is not the whole genomes' worth of data?)" ) ;
		}
		if( m_sample_statistic_filename != "" && m_sample_filename == "" ) {
			m_warnings.push_back( "You are outputting a sample statistic file, but no input sample file has been supplied.\n"
			"   Statistics will be output but the ID fields will be left blank.") ;
		}
		if( m_gen_output_filename == "" && m_sample_output_filename == "" && m_gen_statistic_filename == "" && m_sample_statistic_filename == "" ) {
			m_warnings.push_back( "You have not specified any output files.  This will produce only console output." ) ;
		}
		if( (m_gen_output_filename != "") &&  (m_snp_filter->number_of_subconditions() == 0) && (m_sample_filter->number_of_subconditions() == 0) && (m_gen_statistic_filename == "") && (m_sample_statistic_filename == "") && (m_sample_output_filename == "")) {
			m_warnings.push_back(
			"You have not specified any filters, nor any statistic output files,\n"
			"not a sample output filename.  This will just copy input GEN files\n"
			"to output GEN files (taking into account any file format changes).\n"
			"You can do this faster with gen-convert." ) ;
		}
		if( strings_are_nonempty_and_equal( m_sample_output_filename, m_sample_filename )) {
			m_warnings.push_back( "The input sample file \"" + m_sample_filename + "\" will be overwritten with the output sample file.\n"
				"  (A backup will be taken, but if overwriting is not desired, please use the -os option to choose a different filename.)" ) ;
		}
	}
	
	bool strings_are_nonempty_and_equal( std::string const& left, std::string const& right ) {
		return (!left.empty()) && (!right.empty()) && (left == right) ;
	}
	
	void unsafe_process() {
		process_gen_rows() ;
		process_sample_rows() ;
	}

	void process_gen_rows() {
		backup_file_if_necessary( m_gen_output_filename ) ;
		backup_file_if_necessary( m_gen_statistic_filename ) ;
		m_cout << "Processing SNPs...\n" ;
		Timer timer ;
		open_gen_row_sink() ;
		open_gen_stats_file() ;

		if( genStatisticOutputFile.get() ) { 
			*genStatisticOutputFile << "        " ;
			m_row_statistics.format_column_headers( *genStatisticOutputFile ) << "\n";
		}

		InternalStorageGenRow row ;
		double last_time = -1 ;
		std::size_t number_of_filtered_in_snps = 0 ;
		std::size_t number_of_snps_processed = 0 ;
		for( row.set_number_of_samples( m_number_of_samples_from_gen_file ) ; read_gen_row( row ); row.set_number_of_samples( m_number_of_samples_from_gen_file )) {
			preprocess_gen_row( row ) ;
			process_gen_row( row, ++number_of_snps_processed, number_of_filtered_in_snps ) ;
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
		std::cerr << "Processed " << m_total_number_of_snps << " SNPs in "
			<< std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
		if( m_snp_filter->number_of_subconditions() > 0 ) {
			std::cerr << "(" << number_of_filtered_in_snps << " of " << m_total_number_of_snps << " SNPs passed the filter.)\n" ;
		}
		m_cout << "Post-processing..." << std::flush ;
		timer.restart() ;
		// Close the output gen file(s) now
		close_gen_row_sink() ;
		std::cerr << " (" << std::fixed << std::setprecision(1) << timer.elapsed() << "s)\n\n" ;
	}

	bool read_gen_row( GenRow& row ) {
		return (*gen_row_source) >> row ;
	}

	void preprocess_gen_row( InternalStorageGenRow& row ) const {
		check_gen_row( row ) ;
		row.filter_out_samples_with_indices( m_indices_of_filtered_out_samples ) ;
	}
	
	void check_gen_row( GenRow& row ) const {
		check_gen_row_has_correct_number_of_samples( row ) ;
	}

	void check_gen_row_has_correct_number_of_samples( GenRow& row ) const {
		if( have_sample_file() ) {
			if( row.number_of_samples() != m_number_of_sample_file_rows ) {
				throw GenAndSampleFileMismatchException( "GEN file and sample file have mismatching number of samples." ) ;
			}
		}
	}

	void process_sample_rows() {
		m_cout << "Processing samples...\n" ;
		Timer timer ;

		apply_sample_filter() ;
		assert( m_sample_rows.size() == m_per_column_amounts.size() ) ;

		backup_file_if_necessary( m_sample_output_filename ) ;
		backup_file_if_necessary( m_sample_statistic_filename ) ;
		open_sample_row_sink() ;
		open_sample_stats_file() ;

		if( sampleStatisticOutputFile.get() ) {
			*sampleStatisticOutputFile << "        " ;
			m_sample_statistics.format_column_headers( *sampleStatisticOutputFile ) << "\n" ;
		}

		for( std::size_t i = 0 ; i < m_per_column_amounts.size(); ++i ) {
			SampleRow& sample_row = m_sample_rows[i] ;
			m_sample_statistics.process( sample_row, m_per_column_amounts[i], m_total_number_of_snps ) ;
			if( sampleStatisticOutputFile.get() ) {
				*sampleStatisticOutputFile << std::setw(8) << (i+1) << m_sample_statistics << "\n" ;
			}
			m_sample_statistics.add_to_sample_row( sample_row, "missing" ) ;
			m_sample_statistics.add_to_sample_row( sample_row, "heterozygosity" ) ;
			(*m_sample_row_sink) << sample_row ;
		}

		m_cout << "Processed " << m_number_of_samples_from_gen_file << " samples in " << std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
		if( m_sample_filter->number_of_subconditions() > 0 ) {
			m_cout << "(" << m_sample_rows.size() << " of " << m_number_of_samples_from_gen_file << " samples passed the filter.)\n" ;
		}
	}

	std::string backup_file_if_necessary( std::string const& filename ) {
		try {
			if( filename != "" && exists( filename )) {
				std::string backup_filename = tmpnam(0) ;
				m_cout << "Backing up existing file \"" << filename << "\" as \"" << backup_filename << "\"...\n" ;
				rename( filename, backup_filename ) ;
				return backup_filename ;
			}
		}
		catch( FileException const& e ) {
			m_cout << "A problem occured backing up the file \"" << filename << "\".\n"
				<< e.what() << ".\n" ;
			throw ;
		}
		return "" ;
	}

	void apply_sample_filter() {
		for( std::size_t pre_filter_i = 0, post_filter_i = 0; post_filter_i < m_sample_rows.size(); ++pre_filter_i ) {
			if( sample_row_is_filtered_out( pre_filter_i ) ) {
				if( m_diagnose_sample_filter ) {
						print_sample_filter_diagnostics( m_sample_rows[ post_filter_i ], pre_filter_i ) ;
				}
				m_sample_rows.erase( m_sample_rows.begin() + post_filter_i ) ;
			}
			else {
				++post_filter_i ;
			}
		}
	}

	bool sample_row_is_filtered_out( std::size_t const sample_row_index ) {
		return std::binary_search( m_indices_of_filtered_out_samples.begin(), m_indices_of_filtered_out_samples.end(), sample_row_index ) ;
	}

	void print_sample_filter_diagnostics( SampleRow const& sample_row, std::size_t const sample_row_index ) {
		m_cout
			<< "Filtered out sample row " << sample_row_index << " (" << sample_row.ID1() << " " << sample_row.ID2() << ")"
			<< " because it does not satisfy " ;
		for( std::size_t i = 0, failed_condition_count = 0; i < m_sample_filter->number_of_subconditions() ; ++i ) {
			if( !m_sample_filter->subcondition(i).check_if_satisfied( sample_row )) {
				if( failed_condition_count > 0 ) {
					m_cout << " or " ;
				}
				m_cout << "\"" << m_sample_filter->subcondition(i) << "\"" ;
				++failed_condition_count ;
			}
		}
		m_cout << ".\n" ;
	}

	void process_gen_row( GenRow const& row, std::size_t row_number, std::size_t& number_of_filtered_in_snps ) {
		m_row_statistics.process( row ) ;
		if( m_snp_filter->check_if_satisfied( m_row_statistics )) {
			++number_of_filtered_in_snps ;
			output_gen_row( row ) ;
			output_gen_row_stats( m_row_statistics, row_number ) ;
		}
		else if( m_diagnose_snp_filter ) {
			print_snp_filter_diagnostics( m_row_statistics, row_number ) ;
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

	void print_snp_filter_diagnostics( GenRowStatistics const& row_statistics, std::size_t const row_index ) {
		m_cout
			<< "\rFiltered out snp #" << row_index << " (" << row_statistics.row().SNPID() << " " << row_statistics.row().RSID() << " " << row_statistics.row().SNP_position() << ")"
			<< " because it does not satisfy " ;
		for( std::size_t i = 0, failed_condition_count = 0; i < m_snp_filter->number_of_subconditions() ; ++i ) {
			if( !m_snp_filter->subcondition(i).check_if_satisfied( row_statistics )) {
				if( failed_condition_count > 0 ) {
					m_cout << " or " ;
				}
				m_cout << "\"" << m_snp_filter->subcondition(i) << "\"" ;
				++failed_condition_count ;
			}
		}
		m_cout << ".\n" ;
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
	
	bool have_sample_file() const {
		return m_sample_filename != "" ;
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
	
	std::vector< GenotypeProportions > m_per_column_amounts ;

	std::vector< std::size_t > m_indices_of_filtered_out_samples ;
	
	std::vector< std::string > m_gen_filenames ;
	std::string m_sample_filename ;
	std::string m_gen_output_filename ;
	std::string m_sample_output_filename ;
	std::string m_gen_statistic_filename ;
	std::string m_sample_statistic_filename ;
	
	std::vector< std::string > m_warnings ;
	std::vector< std::string > m_errors ;
	
	std::vector< SampleRow > m_sample_rows ;
	
	bool m_ignore_warnings ;
	bool m_diagnose_sample_filter ;
	bool m_diagnose_snp_filter ;
} ;

int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		GenSelectProcessor::declare_options( options ) ;
		options.process( argc, argv ) ;
    }
	catch( OptionProcessorMutuallyExclusiveOptionsSuppliedException const& e ) {
		std::cerr << "Options \"" << e.first_option()
			<< "\" and \"" << e.second_option()
			<< "\" cannot be supplied at the same time.\n"
			<< "Please consult the documentation, or use \"qc-tool -help\" for more information.\n" ;
		return -1 ;
	}
	catch( OptionProcessorHelpRequestedException const& ) {
	    std::cerr << "Usage: qc-tool <options>\n"
			<< "\nOPTIONS:\n"
            << options
            << "\n" ;
		return 0 ;
	}
    catch( std::exception const& exception ) {
        std::cerr << "!! Error: " << exception.what() << ".\n";
		std::cerr << "Please use \"qc-tool -help\" to see a usage description.\n" ;
        return -1 ;
    }

	// main section

	try	{
		GenSelectProcessor processor( options ) ;
		processor.preprocess() ;
		processor.write_preamble( std::cout ) ;
		processor.process() ;
	}
	catch( ProblemsWereEncountered const& ) {
		return -1 ;
	}
	catch( std::exception const& e )
	{
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}


bool parse_number_from_1_to_100( std::string const& a_string, int* number ) {
	int a_number = parse_integer_in_half_open_range( a_string, 1, 101 ) ;
	if( a_number != 101 ) {
		*number = a_number ;
		std::cout << "number is " << a_number << ".\n" ;
		return true ;
	}
	else {
		return false ;
	}
}

std::vector< std::string > expand_filename_wildcards( std::string const& option_name, std::vector< std::string > const& filenames ) {
	std::vector< std::string > result ;
	for( std::size_t i = 0; i < filenames.size(); ++i ) {
		if( filenames[i].find( '#' ) != std::string::npos ) {
			std::map< int, std::string > results_by_number ;
			std::pair< std::vector< std::string >, std::vector< std::string > > expanded_filename = find_files_matching_path_with_wildcard( filenames[i], '#' ) ;
			if( expanded_filename.first.empty() ) {
				throw OptionValueInvalidException( option_name, filenames, "No file can be found matching filename \"" + filenames[i] + "\"." ) ;
			}

			// We only allow matches corresponding to numbers 1 to 100
			for( std::size_t j = 0; j < expanded_filename.first.size(); ++j ) {
				int number ;
				if( parse_number_from_1_to_100( expanded_filename.second[j], &number )) {
					results_by_number[ number ] = expanded_filename.first[j] ;
				}
			}

			// Copy filenames over (in numerical order)
			for(
				std::map< int, std::string >::const_iterator result_i = results_by_number.begin() ;
				result_i != results_by_number.end() ;
				++result_i
			) {
				result.push_back( result_i->second ) ;
			}
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
