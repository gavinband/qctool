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
#include <boost/bind.hpp>

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
#include "ObjectSource.hpp"
#include "SimpleFileObjectSource.hpp"
#include "SimpleFileObjectSink.hpp"
#include "SNPDataSourceChain.hpp"
#include "SNPDataSinkChain.hpp"
#include "TrivialSNPDataSink.hpp"
#include "TrivialSNPDataSink.hpp"
#include "SampleInputFile.hpp"
#include "SampleOutputFile.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "wildcard.hpp"
#include "string_utils.hpp"
#include "parse_utils.hpp"
#include "progress_bar.hpp"
#include "FileBackupCreator.hpp"
#include "InputToOutputFilenameMapper.hpp"

void check_files_are_readable( std::string const& option_name, std::vector< std::string > const& filenames ) ;
void check_condition_spec( std::string const& option_name, std::string const& condition_spec ) ;

struct QCToolProcessorException: public std::exception
{
	char const* what() const throw() {return "QCToolProcessorException" ; }
} ;

struct GenAndSampleFileMismatchException: public QCToolProcessorException
{
	char const* what() const throw() {return "GenAndSampleFileMismatchException" ; }
} ;

struct NumberOfSamplesMismatchException: public QCToolProcessorException
{
	char const* what() const throw() {return "NumberOfSamplesMismatchException" ; }
} ;

struct SampleFileHasTooFewSamplesException: public QCToolProcessorException
{
	char const* what() const throw() {return "SampleFileHasTooFewSamplesException" ; }
} ;

struct ProblemsWereEncountered: public QCToolProcessorException
{
	char const* what() const throw() {return "ProblemsWereEncountered" ; }
} ;

// Thrown to indicate that the numbers of input and output files specified on the command line differ.
struct QCToolFileCountMismatchError: public QCToolProcessorException
{
	char const* what() const throw() {return "QCToolFileCountMismatchError" ; }
} ;

struct NoMatchingFileFound: public QCToolProcessorException
{
	NoMatchingFileFound( std::string const& filename ): m_filename( filename ) {}
	~NoMatchingFileFound() throw() {}
	char const* what() const throw() { return "NoMatchingFileFound" ; }
	std::string const& filename() const { return m_filename ; }
private:
	std::string m_filename ;
} ;

class OstreamTee: public std::ostream  {
public:

	void add_stream( std::string const& name, std::ostream& stream ) {
		m_streams[ name ] = &stream ;
	}

	template< typename T >
	friend OstreamTee& operator<<( OstreamTee & ostream_tee, T const& t ) ;

	std::ostream& operator[]( std::string const& name ) { return *(m_streams.find( name )->second) ; }

private:
	std::map< std::string, std::ostream* > m_streams ;
} ;

template< typename T >
OstreamTee& operator<<( OstreamTee& ostream_tee, T const& t ) {
	for(
		std::map< std::string, std::ostream* >::const_iterator i = ostream_tee.m_streams.begin() ;
		i != ostream_tee.m_streams.end() ;
		++i
	) {
		(*(i->second)) << t ;
	}

	return ostream_tee ;
}


struct QCToolProcessor: public FileBackupCreator
{
public:
	static void declare_options( OptionProcessor & options ) {

		// Meta-options
		options.set_help_option( "-help" ) ;
		
		// File options	
		options.declare_group( "Data file options" ) ;
	    options[ "-g" ]
	        .set_description( 	"Path of gen file(s) to input.  "
								"To specify several files, either repeat this option or use the numerical wildcard character '#', which "
								"matches numbers from 1 to 100.  For example, \"qc-tool -g myfile_#.gen\" will find all files of "
								"the form \"myfile_N.gen\", where N can be 1, 002, 099, etc." )
	        .set_is_required()
			.set_takes_values()
			.set_maximum_number_of_repeats(23) ;

	    options[ "-og" ]
	        .set_description( 	"Path of gen file(s) to output.  If this option is used, it must appear the same number of times as the -g option. "
	 							"If the corresponding occurence of -g uses a '#' wildcard character, the '#' character can "
								"also be used here to specify numbered output files corresponding to the input files." )
	        .set_takes_values()
			.set_maximum_number_of_repeats(23) ;

	    options[ "-s" ]
	        .set_description( "Path of sample file to input" )
	        .set_takes_single_value()
	        .add_value_checker( &check_files_are_readable ) ;

		options[ "-os" ]
	        .set_description( "Path of sample file to output" )
	        .set_takes_single_value() ;

		// Statistic file options
		options.declare_group( "Statistic file options" ) ;
	    options[ "-snp-stats" ]
	        .set_description( "Output snp-wise statistics to the given file.  If used, it must appear as many times as the -g option.  "
	 							"If the corresponding occurence of -g uses a '#' wildcard character, the '#' character can "
								"also be used here to specify numbered output files corresponding to the input files." )
	        .set_takes_values()
			.set_maximum_number_of_repeats(23) ;

	    options[ "-sample-stats" ]
	        .set_description( "Output sample-wise statistics to the given file." )
	        .set_takes_single_value() ;

		options[ "-snp-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in genstat file.  "
	 						"By default, the columns in this file are: "
							"SNPID, RSID, position, minor_allele, major_allele, MAF, HWE, and missing" )
			.set_takes_single_value()
			.set_default_value( "SNPID, RSID, position, minor_allele, major_allele, MAF, HWE, missing" ) ;

		options[ "-sample-statistics" ]
	        .set_description( "Comma-seperated list of statistics to calculate in samplestat file."
	 						 "  By default, the columns in this file are: ID1, ID2, missing, and heterozygosity.")
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
		options [ "-log" ]
			.set_description( "Path of log file written by qc-tool.")
			.set_takes_single_value()
			.set_default_value( std::string( "qc-tool.log" )) ;
		options [ "-plot" ]
			.set_description( "Path of file to produce plots in.")
			.set_takes_single_value() ;
			
		options.option_excludes_group( "-snp-stats", "SNP filtering options" ) ;
		options.option_excludes_group( "-sample-stats", "Sample filtering options" ) ;
	}

	QCToolProcessor( OptionProcessor const& options )
		: m_options( options )
	{
		m_cout.add_stream( "screen", std::cout ) ;
		write_start_banner( m_cout ) ;
		setup() ;
		preprocess() ;
	}

	~QCToolProcessor()
	{
		write_end_banner( m_cout ) ;
	}
	
private:
	void setup() {
		get_required_filenames() ;
		open_log_file() ;
		try {
			open_m_gen_row_source() ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_cout << "The GEN files specified did not all have the same sample size.\n" ;
			throw ;
		}
		construct_snp_statistics() ;
		construct_sample_statistics() ;
		construct_snp_filter() ;
		construct_sample_filter() ;
		process_other_options() ;
	}

	void get_required_filenames() {
		get_gen_filenames() ;

		if( m_options.check_if_option_was_supplied( "-s" ) ) {
			m_sample_filename = m_options.get_value< std::string >( "-s" ) ;
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
		if( m_options.check_if_option_was_supplied( "-plot" )) {
			m_plot_filename = m_options.get_value< std::string >( "-plot" ) ;
		}
	}

	void get_gen_filenames() {
		assert( m_options.check_if_option_was_supplied( "-g" )) ;
		std::vector< std::string >
			input_gen_filenames_supplied = m_options.get_values< std::string >( "-g" ),
			output_gen_filenames_supplied( input_gen_filenames_supplied.size(), "" ),
			output_snp_stats_filenames_supplied( input_gen_filenames_supplied.size(), "" ) ;

		if( m_options.check_if_option_was_supplied( "-og" ) ) {
			output_gen_filenames_supplied = m_options.get_values< std::string >( "-og" ) ;
		}

		if( m_options.check_if_option_was_supplied( "-snp-stats" ) ) {
			output_snp_stats_filenames_supplied = m_options.get_values< std::string >( "-snp-stats" ) ;
		}

		if( output_gen_filenames_supplied.size() != input_gen_filenames_supplied.size() ) {
			throw QCToolFileCountMismatchError() ;
		}

		if( output_snp_stats_filenames_supplied.size() != input_gen_filenames_supplied.size() ) {
			throw QCToolFileCountMismatchError() ;
		}

		for( std::size_t i = 0; i < input_gen_filenames_supplied.size(); ++i ) {
			m_gen_file_mapper.add_filename_pair( input_gen_filenames_supplied[i], output_gen_filenames_supplied[i] ) ;
			m_snp_stats_file_mapper.add_filename_pair( input_gen_filenames_supplied[i], output_snp_stats_filenames_supplied[i] ) ;
		}
	}

	void open_log_file() {
		m_log_filename = m_options.get_value< std::string >( "-log" ) ;
		backup_file_if_necessary( m_log_filename ) ;
		m_log.open( m_log_filename.c_str() ) ;
		if( m_log.is_open() ) {
			// Make all output go to the log as well.
			m_cout.add_stream( "log", m_log ) ;
		}
	}
	
	void open_m_gen_row_source() {
		Timer timer ;
		
		std::auto_ptr< genfile::SNPDataSourceChain > chain( new genfile::SNPDataSourceChain() ) ;
		chain->set_moved_to_next_source_callback( boost::bind( &QCToolProcessor::move_to_next_output_file, this, _1 )) ;

		for( std::size_t i = 0; i < m_gen_file_mapper.input_files().size(); ++i ) {
			Timer file_timer ;
			m_cout << "(Opening gen file \"" << m_gen_file_mapper.input_files()[i] << "\"...)" << std::flush ;
			try {
				chain->add_source( genfile::SNPDataSource::create( m_gen_file_mapper.input_files()[i] ) ) ;
			}
			catch ( genfile::FileHasTwoConsecutiveNewlinesError const& e ) {
				std::cerr << "\n!!ERROR: a GEN file was specified having two consecutive newlines.\n"
					<< "!! NOTE: popular editors, such as vim and nano, automatically add an extra newline to the file (which you can't see).\n"
					<< "!!     : Please check that each SNP in the file is terminated by a single newline.\n" ;
				throw ;
			}
			m_cout << " (" << file_timer.elapsed() << "s)\n" ;
		}

		m_gen_row_source = chain ;

		if( timer.elapsed() > 1.0 ) {
			m_cout << "Opened " << m_gen_file_mapper.input_files().size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
		}
	}

	void move_to_next_output_file( std::size_t index ) {
		if( index < m_gen_file_mapper.input_files().size() ) {
			if( m_gen_file_mapper.output_filenames().size() > 0 ) {
				if( m_gen_file_mapper.filename_corresponding_to( index ) != m_gen_row_sink->index_of_current_sink() ) {
					m_gen_row_sink->move_to_next_sink() ;
				}
			}
			
			if( m_snp_stats_file_mapper.output_filenames().size() > 0 ) {
				if( m_snp_stats_file_mapper.filename_corresponding_to( index ) != m_current_snp_stats_filename_index ) {
					open_gen_stats_file( ++m_current_snp_stats_filename_index, m_row_statistics ) ;
				}
			}
		}
	}

	void open_gen_row_sink() {
		reset_gen_row_sink() ;
		if( m_gen_file_mapper.output_filenames().size() == 0 ) {
			m_gen_row_sink->add_sink( std::auto_ptr< genfile::SNPDataSink >( new genfile::TrivialSNPDataSink() )) ;
		}
		else {
			for( std::size_t i = 0; i < m_gen_file_mapper.output_filenames().size(); ++i ) {
				m_gen_row_sink->add_sink( genfile::SNPDataSink::create( m_gen_file_mapper.output_filenames()[i] )) ;
			}
		}
	}

	void reset_gen_row_sink() {
		m_gen_row_sink = std::auto_ptr< genfile::SNPDataSinkChain >( new genfile::SNPDataSinkChain() ) ;
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
	
	void reset_sample_statistics_file() {
		sampleStatisticOutputFile.reset() ;
	}

	void reset_gen_statistics_file() {
		genStatisticOutputFile.reset() ;
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
	}
	
public:
	
	void write_start_banner( std::ostream& ) {
		m_cout << "\nWelcome to qc-tool\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}

	void write_end_banner( std::ostream& ) {
		m_cout << "\nThank you for using qc-tool.\n" ;
	}
		
	void write_preamble() {
		m_cout << std::string( 72, '=' ) << "\n\n" ;
		try {
			m_cout << std::setw(30) << "Input SAMPLE file:"
				<< "  \"" << m_sample_filename << "\".\n" ;
			m_cout << std::setw(30) << "Output SAMPLE file:"
				<< "  \"" << m_sample_output_filename << "\".\n" ;
			m_cout << std::setw(30) << "Sample statistic output file:"
				<< "  \"" << m_sample_statistic_filename << "\".\n" ;
			m_cout << "\n" ;

			m_cout << std::setw(30) << "Input GEN file(s):" ;
			for( std::size_t i = 0; i < m_gen_file_mapper.input_files().size(); ++i ) {
				if( i > 0 ) {
					m_cout << std::string( 30, ' ' ) ;
				}
				m_cout << "  (" << std::setw(6) << m_gen_row_source->number_of_snps_in_source( i ) << " snps)  " ;
				m_cout << "\"" << m_gen_file_mapper.input_files()[i] << "\"" ;				
				m_cout << "\n" ;
			}
			m_cout << std::string( 30, ' ' ) << "  (total " << m_gen_row_source->total_number_of_snps() << " snps in input).\n" ;
			if( m_gen_file_mapper.input_files().size() > 1 ) {
				m_cout << "\n" ;
			}

			m_cout << std::setw(30) << "Output GEN file(s):" ;
			if( m_gen_file_mapper.output_filenames().empty() ) {
				m_cout << "  (n/a)\n" ;
			}
			else {
				for( std::size_t i = 0; i < m_gen_file_mapper.output_filenames().size(); ++i ) {
					if( i > 0 ) {
						m_cout << "\n" << std::string( 30, ' ' ) ;
					}
					m_cout << "  \"" << m_gen_file_mapper.output_filenames()[i] << "\"" ;				
				}
				m_cout << "\n" ;
			}

			m_cout << std::setw(30) << "SNP statistic output file(s):" ;
			for( std::size_t i = 0; i < m_snp_stats_file_mapper.output_filenames().size(); ++i ) {
				if( i > 0 ) {
					m_cout << std::string( 30, ' ' ) ;
				}
				m_cout << "  \"" << m_snp_stats_file_mapper.output_filenames()[i] << "\"\n" ;
			}
			m_cout << "\n" ;
			m_cout << std::setw(30) << "Sample filter:" 
				<< "  " << *m_sample_filter << ".\n" ;
			m_cout << std::setw(30) << "SNP filter:"
				<< "  " << *m_snp_filter << ".\n" ;
			m_cout << "\n" ;

			m_cout << std::setw(30) << "# of samples in input files:"
				<< "  " << m_gen_row_source->number_of_samples() << ".\n" ;
			m_cout << std::setw(30) << "# of samples after filtering:"
				<< "  " << m_gen_row_source->number_of_samples() - m_indices_of_filtered_out_samples.size()
				<< " (" << m_indices_of_filtered_out_samples.size()
				<< " filtered out).\n" ;
			
			m_cout << "\n" << std::string( 72, '=' ) << "\n\n" ;

			if( !m_errors.empty() ) {
				for( std::size_t i = 0; i < m_errors.size(); ++i ) {
					m_cout << "!! ERROR: " << m_errors[i] << "\n\n" ;
				}
				m_cout << "!! Please correct the above errors and re-run qc-tool.\n" ;
				throw ProblemsWereEncountered() ;
			}


			if( !m_warnings.empty() ) {
				for( std::size_t i = 0; i < m_warnings.size(); ++i ) {
					m_cout << "!! WARNING: " << m_warnings[i] << "\n\n" ;
				}
				if( m_ignore_warnings ) {
					m_cout << "!! Warnings were encountered, but proceeding anyway as -force was supplied.\n" ;
					m_cout << "\n" << std::string( 72, '=' ) << "\n\n" ;
				}
				else {
					m_cout << "!! Warnings were encountered.  To proceed anyway, please run again with the -force option.\n" ;
					throw ProblemsWereEncountered() ;
				}
			}
		}
		catch (...) {
			throw ;
		}
	}

	void write_postamble() {
		m_cout << std::string( 72, '=' ) << "\n\n" ;
		try {
			if( backed_up_files().size() > 0 ) {
				m_cout << std::setw(36) << "I took backups of the following files before overwriting:\n" ;
				std::size_t max_length = 0u ;
				for(
					std::map< std::string, std::string >::const_iterator i = backed_up_files().begin() ;
					i != backed_up_files().end() ;
					++i
				) {
					max_length = std::max( max_length, i->first.size() ) ;
				}

				for(
					std::map< std::string, std::string >::const_iterator i = backed_up_files().begin() ;
					i != backed_up_files().end() ;
					++i
				) {
					m_cout << "  " << std::setw( max_length + 2 ) << std::left << ("\"" + i->first + "\"") << " to \"" << i->second << "\"\n" ;
				}

				m_cout << "\n" ;
				m_cout << std::string( 72, '=' ) << "\n\n" ;
			}

			m_cout << std::setw(36) << "Number of SNPs in input file(s):"
				<< "  " << m_gen_row_source->total_number_of_snps() << ".\n" ;
			if( m_snp_filter->number_of_subconditions() > 0 ) {
				for(
					std::map< std::size_t, std::size_t >::const_iterator i = m_snp_filter_failure_counts.begin() ;
					i != m_snp_filter_failure_counts.end(); 
					++i
				) {
					m_cout << std::setw(36) << ("...which failed \"" + to_string( m_snp_filter->subcondition( i->first )) + "\":")
						<< "  " << i->second << ".\n" ;
				}

				m_cout << std::setw(36) << "(total failures:"
					<< "  " << m_gen_row_source->total_number_of_snps() - m_number_of_filtered_in_snps << ").\n" ;
			}
			
			m_cout << "\n" ;

			m_cout << std::setw(36) << "Number of samples in input file(s):"
				<< "  " << m_gen_row_source->number_of_samples() << ".\n" ;
			if( m_sample_filter->number_of_subconditions() > 0 ) {
				for(
					std::map< std::size_t, std::size_t >::const_iterator i = m_sample_filter_failure_counts.begin() ;
					i != m_sample_filter_failure_counts.end(); 
					++i
				) {
					m_cout << std::setw(36) << ("...which failed \"" + to_string( m_sample_filter->subcondition( i->first )) + "\":")
						<< "  " << i->second << ".\n" ;
				}
				m_cout << std::setw(36) << "(total failures:"
					<< "  " << m_gen_row_source->number_of_samples() - m_sample_rows.size() << ").\n" ;
			}

			m_cout << "\n" ;

			if( m_gen_file_mapper.output_filenames().size() > 0 ) {
				m_cout << std::setw(30) << "Output GEN files:" ;
				for( std::size_t i = 0; i < m_gen_file_mapper.output_filenames().size(); ++i ) {
					if( i > 0 ) {
						m_cout << std::string( 30, ' ' ) ;
					}
					if( m_gen_row_sink.get() ) {
						m_cout << "  (" << std::setw(6) << m_gen_row_sink->sink(i).number_of_snps_written() << " snps)  " ;
					}
					m_cout << "\"" << m_gen_file_mapper.input_files()[i] << "\"\n" ;
				}
				m_cout << std::string( 30, ' ' ) << "  (total " << m_gen_row_sink->number_of_snps_written() << " snps).\n" ;
			}
			if( m_sample_output_filename != "" ) {
				m_cout << std::setw(36) << "Output SAMPLE files:"
					<< "  \"" << m_sample_output_filename << "\""
					<< "  (" << m_sample_rows.size() << " samples)\n" ;
			}
			if( m_snp_stats_file_mapper.output_filenames().size() > 0 ) {
				m_cout << std::setw(30) << "SNP statistic output file(s):" ;
				for( std::size_t i = 0; i < m_snp_stats_file_mapper.output_filenames().size(); ++i ) {
					if( i > 0 ) {
						m_cout << std::string( 30, ' ' ) ;
					}
					m_cout << "  \"" << m_snp_stats_file_mapper.output_filenames()[i] << "\"\n" ;
				}
			}
			if( m_sample_statistic_filename != "" ) {
				m_cout << std::setw(36) << "Sample statistic output file:"
					<< "  \"" << m_sample_statistic_filename << "\".\n" ;
			}

			m_cout[ "screen" ] << std::setw( 36 ) << "\nMore details are in the log file:"
				<< "  \"" << m_log_filename << "\".\n" ;
			m_cout << std::string( 72, '=' ) << "\n\n" ;
		}
		catch (...) {
			throw ;
		}
	}

	void preprocess() {
		open_sample_row_source() ;
		preprocess_sample_rows() ;
		check_for_errors_and_warnings() ;
	}

	void process() {
		try {
			write_preamble() ;
			unsafe_process() ;
			write_postamble() ;
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
			if( m_number_of_sample_file_rows != m_gen_row_source->number_of_samples() ) {
				throw NumberOfSamplesMismatchException() ;
			}
		}
		else {
			m_sample_rows.resize( m_gen_row_source->number_of_samples() ) ;
		}
	}

	void check_for_errors_and_warnings() {
		check_for_errors() ;
		check_for_warnings() ;
	}

	void check_for_errors() {
		if( m_gen_file_mapper.input_files().size() == 0 ) {
			m_errors.push_back( "At least one GEN input file must be supplied." ) ;
		}

		for( std::size_t i = 0; i < m_gen_file_mapper.input_files().size(); ++i ) {
			for( std::size_t j = 0; j < m_gen_file_mapper.output_filenames().size(); ++j ) {
				if( strings_are_nonempty_and_equal( m_gen_file_mapper.output_filenames()[j], m_gen_file_mapper.input_files()[i] )) {
					m_errors.push_back( "Output GEN file \"" + m_gen_file_mapper.output_filenames()[j] +"\" also specified as input GEN file." ) ;
					break ;
				}
			}
			if( strings_are_nonempty_and_equal( m_sample_output_filename, m_gen_file_mapper.input_files()[i] )) {
				m_errors.push_back( "Output SAMPLE file \"" + m_sample_output_filename +"\" also specified as input GEN file." ) ;
				break ;
			}
			for( std::size_t j = 0; j < m_snp_stats_file_mapper.output_filenames().size(); ++j ) {
				if( strings_are_nonempty_and_equal( m_snp_stats_file_mapper.output_filenames()[j], m_gen_file_mapper.input_files()[i] )) {
					m_errors.push_back( "Output GEN statistic file \"" + m_snp_stats_file_mapper.output_filenames()[j] +"\" also specified as input GEN file." ) ;
					break ;
				}
			}
			if( strings_are_nonempty_and_equal( m_sample_statistic_filename, m_gen_file_mapper.input_files()[i] )) {
				m_errors.push_back( "Output SAMPLE statistic file \"" + m_sample_statistic_filename +"\" also specified as input GEN file." ) ;
				break ;
			}
		}
		for( std::size_t j = 0; j < m_gen_file_mapper.output_filenames().size(); ++j ) {
			if( strings_are_nonempty_and_equal( m_sample_output_filename, m_gen_file_mapper.output_filenames()[j] )) {
				m_errors.push_back( "The GEN and SAMPLE output filenames must differ." ) ;
			}
		}
		for( std::size_t j = 0; j < m_snp_stats_file_mapper.output_filenames().size(); ++j ) {
			if( strings_are_nonempty_and_equal( m_snp_stats_file_mapper.output_filenames()[j], m_sample_statistic_filename )) {
				m_errors.push_back( "The gen statistic and sample statistic filenames must differ." ) ;
			}
		}
		if( m_sample_filename == "" && m_sample_filter->number_of_subconditions() != 0 ) {
			m_errors.push_back( "To filter on samples, please supply an input sample file." ) ;
		}
	}
	
	void check_for_warnings() {
		if( (m_sample_output_filename != "" || m_sample_statistic_filename != "") && m_gen_file_mapper.input_files().size() != 23 ) {
			m_warnings.push_back( "You are outputting a sample or sample statistic file, but the number of gen files is not 23.\n"
			"   (I suspect there is not the whole genomes' worth of data?)" ) ;
		}
		if( m_sample_statistic_filename != "" && m_sample_filename == "" ) {
			m_warnings.push_back( "You are outputting a sample statistic file, but no input sample file has been supplied.\n"
			"   Statistics will be output but the ID fields will be left blank.") ;
		}
		if( m_gen_file_mapper.output_filenames().size() == 0 && m_sample_output_filename == "" && m_snp_stats_file_mapper.output_filenames().size() == 0 && m_sample_statistic_filename == "" ) {
			m_warnings.push_back( "You have not specified any output files.  This will produce only console output." ) ;
		}
		if( (m_gen_file_mapper.output_filenames().size() > 0) &&  (m_snp_filter->number_of_subconditions() == 0) && (m_sample_filter->number_of_subconditions() == 0)) {
			m_warnings.push_back( "You have specified output GEN files, but no filters.\n"
			 	"The output GEN files will contain the same data as the input ones.") ;
		}
		if(m_gen_file_mapper.output_filenames().size() == 0 && ((m_snp_filter->number_of_subconditions() > 0) || (m_sample_filter->number_of_subconditions() > 0))) {
			m_warnings.push_back( "You have not specified filters, but no output GEN files." ) ;
		}
		if(m_sample_output_filename == "" && (m_sample_filter->number_of_subconditions() > 0)) {
			m_warnings.push_back( "You have not specified sample filters, but no output sample file.\n" ) ;
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
		construct_plots() ;
	}

	void process_gen_rows() {
		for( std::size_t i = 0; i < m_gen_file_mapper.output_filenames().size(); ++i ) {
			backup_file_if_necessary( m_gen_file_mapper.output_filenames()[i] ) ;
		}
		for( std::size_t i = 0; i < m_snp_stats_file_mapper.output_filenames().size(); ++i ) {
			backup_file_if_necessary( m_snp_stats_file_mapper.output_filenames()[i] ) ;
		}
		m_cout << "Processing SNPs...\n" ;
		Timer timer ;
		open_gen_row_sink() ;
		if( m_snp_stats_file_mapper.output_filenames().size() > 0 ) {
			open_gen_stats_file( 0, m_row_statistics ) ;
		}

		InternalStorageGenRow row ;
		double last_time = -1 ;
		m_number_of_filtered_in_snps = 0 ;
		std::size_t number_of_snps_processed = 0 ;
		for( row.set_number_of_samples( m_gen_row_source->number_of_samples() ) ; read_gen_row( row ); row.set_number_of_samples( m_gen_row_source->number_of_samples() )) {
			preprocess_gen_row( row ) ;
			process_gen_row( row, ++number_of_snps_processed, m_number_of_filtered_in_snps ) ;
			accumulate_per_column_amounts( row, m_per_column_amounts ) ;
			double time_now = timer.elapsed() ;
			if( (time_now - last_time > 1.0) || (number_of_snps_processed == m_gen_row_source->total_number_of_snps()) ) {
				m_cout[ "screen" ]
					<< "\r"
					<< get_progress_bar( 30, static_cast< double >( number_of_snps_processed ) / m_gen_row_source->total_number_of_snps() )
					<< " (" << number_of_snps_processed << " / " << m_gen_row_source->total_number_of_snps()
					<< ", " << std::fixed << std::setprecision(1) << time_now << "s)"
					<< std::flush ;
				last_time = time_now ;
			}
		}
		m_cout << "\n" ;
		assert( number_of_snps_processed = m_gen_row_source->total_number_of_snps() ) ;
		m_cout << "Processed " << m_gen_row_source->total_number_of_snps() << " SNPs in "
			<< std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
		if( m_snp_filter->number_of_subconditions() > 0 ) {
			m_cout << "(" << m_number_of_filtered_in_snps << " of " << m_gen_row_source->total_number_of_snps() << " SNPs passed the filter.)\n" ;
		}
		/*
		m_cout << "Post-processing..." << std::flush ;
		timer.restart() ;
		// Close the output gen file(s) now
		reset_gen_row_sink() ;
		m_cout << " (" << std::fixed << std::setprecision(1) << timer.elapsed() << "s)\n\n" ;
		*/
	}

	void open_gen_stats_file( std::size_t index, GenRowStatistics const& row_statistics ) {
		assert( index < m_snp_stats_file_mapper.output_filenames().size()) ;
		m_current_snp_stats_filename_index = index ;
		genStatisticOutputFile = open_file_for_output( m_snp_stats_file_mapper.output_filenames()[ index ] ) ;
		if( genStatisticOutputFile.get() ) { 
			*genStatisticOutputFile << "        " ;
			m_row_statistics.format_column_headers( *genStatisticOutputFile ) << "\n";
		}

	}


	bool read_gen_row( GenRow& row ) {
		return row.read_from_source( *m_gen_row_source ) ;
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
				throw GenAndSampleFileMismatchException() ;
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
			m_sample_statistics.process( sample_row, m_per_column_amounts[i], m_gen_row_source->total_number_of_snps() ) ;
			if( sampleStatisticOutputFile.get() ) {
				*sampleStatisticOutputFile << std::setw(8) << (i+1) << m_sample_statistics << "\n" ;
			}
			m_sample_statistics.add_to_sample_row( sample_row, "missing" ) ;
			m_sample_statistics.add_to_sample_row( sample_row, "heterozygosity" ) ;
			(*m_sample_row_sink) << sample_row ;
		}

		m_cout << "Processed " << m_gen_row_source->number_of_samples() << " samples in " << std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
		if( m_sample_filter->number_of_subconditions() > 0 ) {
			m_cout << "(" << m_sample_rows.size() << " of " << m_gen_row_source->number_of_samples() << " samples passed the filter.)\n" ;
		}
		m_cout << "\n" ;
	}

	void apply_sample_filter() {
		for( std::size_t pre_filter_i = 0, post_filter_i = 0; post_filter_i < m_sample_rows.size(); ++pre_filter_i ) {
			if( sample_row_is_filtered_out( pre_filter_i ) ) {
				do_sample_filter_diagnostics( m_sample_rows[ post_filter_i ], pre_filter_i ) ;
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

	void do_sample_filter_diagnostics( SampleRow const& sample_row, std::size_t const sample_row_index ) {
		std::ostream& log = m_cout["log"] ;
		log
			<< "Filtered out sample row " << sample_row_index << " (" << sample_row.ID1() << " " << sample_row.ID2() << ")"
			<< " because it does not satisfy " ;
		for( std::size_t i = 0, failed_condition_count = 0; i < m_sample_filter->number_of_subconditions() ; ++i ) {
			if( !m_sample_filter->subcondition(i).check_if_satisfied( sample_row )) {
				++m_sample_filter_failure_counts[ i ] ;
				if( failed_condition_count > 0 ) {
					log << " or " ;
				}
				log << "\"" << m_sample_filter->subcondition(i) << "\"" ;
				++failed_condition_count ;
			}
		}
		log << ".\n" ;
	}

	void process_gen_row( GenRow const& row, std::size_t row_number, std::size_t& m_number_of_filtered_in_snps ) {
		m_row_statistics.process( row ) ;
		if( m_snp_filter->check_if_satisfied( m_row_statistics )) {
			++m_number_of_filtered_in_snps ;
			output_gen_row( row ) ;
			output_gen_row_stats( m_row_statistics, row_number ) ;
		}
		else {
			do_snp_filter_diagnostics( m_row_statistics, row_number ) ;
		}
	}

	void output_gen_row( GenRow const& row ) {
		row.write_to_sink( *m_gen_row_sink ) ;
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

	void do_snp_filter_diagnostics( GenRowStatistics const& row_statistics, std::size_t const row_index ) {
		std::ostream& log = m_cout["log"] ;
		log
			<< "Filtered out snp #" << row_index << " (" << row_statistics.row().SNPID() << " " << row_statistics.row().RSID() << " " << row_statistics.row().SNP_position() << ")"
			<< " because it does not satisfy " ;
		for( std::size_t i = 0, failed_condition_count = 0; i < m_snp_filter->number_of_subconditions() ; ++i ) {
			if( !m_snp_filter->subcondition(i).check_if_satisfied( row_statistics )) {

				++m_snp_filter_failure_counts[ i ] ;

				if( failed_condition_count > 0 ) {
					log << " or " ;
				}
				log << "\"" << m_snp_filter->subcondition(i) << "\"" ;
				++failed_condition_count ;
			}
		}
		log << ".\n" ;
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
	
	void construct_plots() {
		// not implemented.
	}
	
private:
	
	std::string m_log_filename ;
	std::ofstream m_log ;
	OstreamTee m_cout ;
	
	std::auto_ptr< genfile::SNPDataSourceChain > m_gen_row_source ;
	std::auto_ptr< genfile::SNPDataSinkChain > m_gen_row_sink ;
	std::auto_ptr< ObjectSource< SampleRow > > m_sample_row_source ;
	std::auto_ptr< ObjectSink< SampleRow > > m_sample_row_sink ;
	OUTPUT_FILE_PTR genStatisticOutputFile ;
	OUTPUT_FILE_PTR sampleStatisticOutputFile ;
	OptionProcessor const& m_options ;
	
	GenRowStatistics m_row_statistics ;
	SampleRowStatistics m_sample_statistics ;
	std::auto_ptr< AndRowCondition > m_snp_filter ;
	std::auto_ptr< AndRowCondition > m_sample_filter ;
	
	std::size_t m_number_of_sample_file_rows ;
	std::size_t m_number_of_filtered_in_snps ;
	
	std::map< std::size_t, std::size_t > m_snp_filter_failure_counts ;
	std::map< std::size_t, std::size_t > m_sample_filter_failure_counts ;
	
	std::vector< GenotypeProportions > m_per_column_amounts ;

	std::vector< std::size_t > m_indices_of_filtered_out_samples ;
		
	InputToOutputFilenameMapper m_gen_file_mapper ;
	InputToOutputFilenameMapper m_snp_stats_file_mapper ;
	std::size_t m_current_snp_stats_filename_index ;

	std::string m_sample_filename ;
	std::string m_sample_output_filename ;
	std::string m_sample_statistic_filename ;
	std::string m_plot_filename ;

	std::vector< std::string > m_warnings ;
	std::vector< std::string > m_errors ;
	
	std::vector< SampleRow > m_sample_rows ;
	
	bool m_ignore_warnings ;
} ;

int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		QCToolProcessor::declare_options( options ) ;
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
		QCToolProcessor processor( options ) ;
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
