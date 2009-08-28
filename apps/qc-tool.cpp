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
#include "CmdLineOptionProcessor.hpp"
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
#include "OstreamTee.hpp"

namespace globals {
	std::string const program_name = "qc-tool" ;
}

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

struct SNPPositionSink: public genfile::SNPDataSink
{
	SNPPositionSink( std::string const& filename ) {
		setup( filename ) ;
	}
	
	operator bool() const { return *m_stream_ptr ; }

	void write_snp_impl(
		uint32_t,
		std::string,
		std::string,
		uint32_t SNP_position,
		char,
		char,
		GenotypeProbabilityGetter const&,
		GenotypeProbabilityGetter const&,
		GenotypeProbabilityGetter const&
	) {
		stream() << SNP_position << "\n" ;
	} ;
	
private:
	
	std::ostream& stream() { return *m_stream_ptr ; }
	
	void setup( std::string const& filename ) {
		m_stream_ptr = open_file_for_output( filename ) ;
	}
	
	OUTPUT_FILE_PTR m_stream_ptr ;
} ;


struct QCToolOptionProcessor: public CmdLineOptionProcessor
{
public:
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( OptionProcessor& options ) {
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
	        .set_description( 	"Override the auto-generated path(s) of the output gen file(s).  "
								"If this option is supplied, it must appear the same number of times as the -g option. "
	 							"If the corresponding occurence of -g uses a '#' wildcard character, the '#' character can "
								"also be used here to specify numbered output files corresponding to the input files." )
	        .set_takes_values()
			.set_maximum_number_of_repeats(23) ;

	    options[ "-s" ]
	        .set_description( "Path of sample file to input" )
			.set_takes_single_value() ;

		options[ "-os" ]
	        .set_description( "Override the auto-generated path of the output sample file.  " )
	        .set_takes_single_value() ;

		// Statistic file options
		options.declare_group( "Statistic calculation options" ) ;
	    options[ "-snp-stats" ]
			.set_description( "Calculated and output snp-wise statistics." ) ;
	    options[ "-snp-stats-file" ]
	        .set_description( 	"Override the auto-generated path(s) of the files in which snp-wise statistics will be output.  "
								"If used, this option must appear as many times as the -g option.  "
	 							"If the corresponding occurence of -g uses a '#' wildcard character, the '#' character can "
								"also be used here to specify numbered output files corresponding to the input files." )
	        .set_takes_values()
			.set_maximum_number_of_repeats(23) ;

		options[ "-snp-stats-columns" ]
	        .set_description( "Comma-seperated list of columns to output in the snp-wise statistics file.  "
	 						"By default, the columns are: "
							"SNPID, RSID, position, minor_allele, major_allele, MAF, HWE, and missing" )
			.set_takes_single_value()
			.set_default_value( "SNPID, RSID, position, minor_allele, major_allele, MAF, HWE, missing" ) ;

	    options[ "-sample-stats" ]
			.set_description( "Calculate and output sample-wise statistics." ) ;
		options[ "-sample-stats-columns" ]
	        .set_description( "Comma-seperated list of statistics to output in the sample-wise statistics file."
	 						 "  By default, the columns are: ID1, ID2, missing, and heterozygosity.")
			.set_takes_single_value()
			.set_default_value( std::string("ID1, ID2, missing, heterozygosity") ) ;
	    options[ "-sample-stats-file" ]
	        .set_description( 	"Override the auto-generated path of the file in which sample-wise statistics will be output." )
	        .set_takes_single_value() ;

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
			.set_description( "Filter out SNPs whose SNP ID or RSID does not lie in the given file(s).")
			.set_takes_values() 
			.set_number_of_values_per_use( 1 ) ;
		options[ "-snp-excl-list"]
			.set_description( "Filter out SNPs whose SNP ID or RSID lies in the given file(s).")
			.set_takes_values() 
			.set_number_of_values_per_use( 1 ) ;
		options [ "-write-excl-list" ]
			.set_description( "Don't apply the filter; instead, write files containing the positions of SNPs that would be filtered out."
			"  These files are suitable for use as the input to -snp-incl-list option on a subsequent run." ) ;
		options [ "-write-excl-list-file" ]
			.set_description( "Don't apply the filter; instead, write files containing the positions of SNPs that would be filtered out."
			"  These files are suitable for use as the input to -snp-incl-list option on a subsequent run." ) ;
		options.option_excludes_option( "-write-excl-list", "-og" ) ;

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

		// Other options
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
	
	void process( int argc, char** argv ) {
		CmdLineOptionProcessor::process( argc, argv ) ;
		process_filenames() ;
	}

	std::string const& input_sample_filename() const { return m_input_sample_filename ; }
	std::string const& output_sample_filename() const { return m_output_sample_filename ; }
	std::string const& output_sample_stats_filename() const { return m_sample_statistic_filename ; }
	std::string const& plot_filename() const { return m_plot_filename ; }
	std::string const& log_filename() const { return m_log_filename ; }
	InputToOutputFilenameMapper gen_filename_mapper() const { return m_gen_file_mapper ; }
	InputToOutputFilenameMapper snp_stats_filename_mapper() const { return m_snp_stats_file_mapper ; }
	InputToOutputFilenameMapper output_snp_excl_filename_mapper() const { return m_output_snp_excl_file_mapper ; }
	std::vector< std::string > row_statistics_specs() const {
		return split_and_strip_discarding_empty_entries( get_value< std::string >( "-snp-stats-columns" ), "," ) ;
	}
	std::vector< std::string > sample_statistics_specs() const {
		return split_and_strip_discarding_empty_entries( get_value< std::string >( "-sample-stats-columns" ), "," ) ;
	}
	
private:
	
	void process_filenames() {
		get_snp_related_filenames() ;
		get_sample_related_filenames() ;
		if( check_if_option_was_supplied( "-plot" )) {
			m_plot_filename = get_value< std::string >( "-plot" ) ;
		}
	}

	void get_snp_related_filenames() {
		assert( check_if_option_was_supplied( "-g" )) ;
		std::vector< std::string >
			input_gen_filenames_supplied = get_values< std::string >( "-g" ),
			output_gen_filenames_supplied = construct_output_gen_filenames( input_gen_filenames_supplied ),
			output_snp_stats_filenames_supplied = construct_snp_stats_filenames( input_gen_filenames_supplied ),
			output_snp_excl_filenames_supplied = construct_output_snp_excl_list_filenames( input_gen_filenames_supplied ) ;

		for( std::size_t i = 0; i < input_gen_filenames_supplied.size(); ++i ) {
			m_gen_file_mapper.add_filename_pair( input_gen_filenames_supplied[i], output_gen_filenames_supplied[i] ) ;
			m_snp_stats_file_mapper.add_filename_pair( input_gen_filenames_supplied[i], output_snp_stats_filenames_supplied[i] ) ;
			m_output_snp_excl_file_mapper.add_filename_pair( input_gen_filenames_supplied[i], output_snp_excl_filenames_supplied[i] ) ;
		}
	}

	void get_sample_related_filenames() {
		if( check_if_option_was_supplied( "-s" ) ) {
			m_input_sample_filename = get_value< std::string >( "-s" ) ;
		}

		if( check_if_option_was_supplied( "-sample-stats" ) ) {
			if( check_if_option_was_supplied( "-sample-stats-file" )) {
				m_sample_statistic_filename = get_value< std::string >( "-sample-stats-file" ) ;
			}
			else {
				m_sample_statistic_filename = strip_sample_file_extension_if_present( m_input_sample_filename ) + ".sample-stats";
			}

			if( check_if_option_was_supplied( "-os" )) {
				m_output_sample_filename = get_value< std::string >( "-os" ) ;
			}
			else {
				m_output_sample_filename = m_input_sample_filename ;
			}
		}

		if( check_if_option_was_supplied_in_group( "Sample filtering options" )) {
			if( check_if_option_was_supplied( "-os" )) {
				m_output_sample_filename = get_value< std::string >( "-os" ) ;
			}
			else {
				m_output_sample_filename = strip_sample_file_extension_if_present( m_input_sample_filename ) + ".fltrd.sample" ;
			}
		}
	}

	std::string strip_sample_file_extension_if_present( std::string filename ) {
		if( filename.size() >= 7 && filename.substr( filename.size() - 7, 7 ) == ".sample" ) {
			filename.resize( filename.size() - 7 ) ;
		}
		return filename ;
	}

	std::vector< std::string > construct_output_gen_filenames( std::vector< std::string > const& input_gen_filenames_supplied ) {
		std::vector< std::string > result( input_gen_filenames_supplied.size(), "" ) ;
		if( check_if_option_was_supplied( "-og" ) ) {
			result = get_values< std::string >( "-og" ) ;
		} else if( check_if_option_was_supplied_in_group( "SNP filtering options" ) || check_if_option_was_supplied_in_group( "Sample filtering options" )) {
			for( std::size_t i = 0; i < input_gen_filenames_supplied.size(); ++i ) {
				result[i]
					= genfile::strip_gen_file_extension_if_present( input_gen_filenames_supplied[i] )
					+ ".fltrd-in"
					+ genfile::get_gen_file_extension_if_present( input_gen_filenames_supplied[i] ) ;
			}
		}
		if( result.size() != input_gen_filenames_supplied.size() ) {
			throw QCToolFileCountMismatchError() ;
		}
		return result ;
	}

	std::vector< std::string > construct_snp_stats_filenames( std::vector< std::string > const& input_gen_filenames_supplied ) {
		std::vector< std::string > result( input_gen_filenames_supplied.size(), "" ) ;
		if( check_if_option_was_supplied( "-snp-stats" ) ) {
			if( check_if_option_was_supplied( "-snp-stats-file" )) {
				result = get_values< std::string >( "-snp-stats-file" ) ;
			}
			else {
				for( std::size_t i = 0; i < input_gen_filenames_supplied.size() ; ++i ) {
					result[i] = genfile::strip_gen_file_extension_if_present( input_gen_filenames_supplied[i] ) + ".snp-stats" ;
				}
			}
		}
		if( result.size() != input_gen_filenames_supplied.size() ) {
			throw QCToolFileCountMismatchError() ;
		}
		return result ;
	}

	std::vector< std::string > construct_output_snp_excl_list_filenames( std::vector< std::string > const& input_gen_filenames_supplied ) {
		std::vector< std::string > result( input_gen_filenames_supplied.size(), "" ) ;
		if( check_if_option_was_supplied( "-write-excl-list" ) ) {
			if( check_if_option_was_supplied( "-write-excl-list-file" )) {
				result = get_values< std::string >( "-write-excl-list-file" ) ;
			}
			else {
				for( std::size_t i = 0; i < input_gen_filenames_supplied.size() ; ++i ) {
					result[i] = genfile::strip_gen_file_extension_if_present( input_gen_filenames_supplied[i] ) + ".snp-excl" ;
				}
			}
		}
		if( result.size() != input_gen_filenames_supplied.size() ) {
			throw QCToolFileCountMismatchError() ;
		}
		return result ;
	}

	
private:	
	std::string m_input_sample_filename ;
	std::string m_output_sample_filename ;
	std::string m_sample_statistic_filename ;
	std::string m_plot_filename ;
	std::string m_log_filename ;
	InputToOutputFilenameMapper m_gen_file_mapper ;
	InputToOutputFilenameMapper m_snp_stats_file_mapper ;
	InputToOutputFilenameMapper m_output_snp_excl_file_mapper ;
} ;


struct QCToolContext
{
	typedef genfile::SNPDataSource SNPDataSource ;
	typedef genfile::SNPDataSink SNPDataSink ;
	
	virtual ~QCToolContext() {}
	
	virtual SNPDataSource& snp_data_source() const = 0 ;
	virtual SNPDataSink& fltrd_in_snp_data_sink() const = 0 ;
	virtual SNPDataSink& fltrd_out_snp_data_sink() const = 0 ;
	virtual ObjectSink< SampleRow >& fltrd_in_sample_sink() const = 0 ;
	virtual ObjectSink< SampleRow >& fltrd_out_sample_sink() const = 0 ;
	virtual std::ostream& snp_stats_sink() const = 0 ;
	virtual std::ostream& sample_stats_sink() const = 0 ;
	virtual AndRowCondition& snp_filter() const = 0 ;
	virtual AndRowCondition& sample_filter() const = 0 ;
	virtual std::vector< std::size_t >& snp_filter_failure_counts() = 0 ;
	virtual std::vector< std::size_t >& sample_filter_failure_counts() = 0 ;
	virtual GenRowStatistics& snp_statistics() = 0 ;
	virtual SampleRowStatistics& sample_statistics() = 0 ;
	virtual std::vector< SampleRow >& sample_rows() = 0 ;

	virtual void print_progress_if_necessary() = 0 ;

	virtual OstreamTee& logger() = 0 ;
} ;


struct QCToolCmdLineContext: public QCToolContext
{
	QCToolCmdLineContext( int argc, char** argv ) {
		timestamp() ;
		write_start_banner() ;
		m_options.process( argc, argv ) ;
		setup() ;
		timestamp() ;
	}
	
	~QCToolCmdLineContext() {
		write_end_banner() ;
	}
	
	SNPDataSource& snp_data_source() const {
		return *m_snp_data_source ;
	}

	SNPDataSink& fltrd_in_snp_data_sink() const {
		return *m_fltrd_in_snp_data_sink ;
	}
	
	SNPDataSink& fltrd_out_snp_data_sink() const {
		return *m_fltrd_out_snp_data_sink ;
	}
	
	std::vector< SampleRow >& sample_rows() {
		return m_sample_rows ;
	}

	ObjectSink< SampleRow >& fltrd_in_sample_sink() const {
		return *m_fltrd_in_sample_sink ;
	}

	ObjectSink< SampleRow >& fltrd_out_sample_sink() const {
		return *m_fltrd_out_sample_sink ;
	}

	std::ostream& snp_stats_sink() const {
		assert( m_snp_stats_file.get() ) ;
		return *m_snp_stats_file ;
	}

	std::ostream& sample_stats_sink() const {
		assert( m_sample_stats_file.get() ) ;
		return *m_sample_stats_file ;
	}
	
	AndRowCondition& snp_filter() const {
		return *m_snp_filter ;
	}
	
	AndRowCondition& sample_filter() const {
		return *m_sample_filter ;
	}
	
	GenRowStatistics& snp_statistics() {
		return m_snp_statistics ;
	}

	SampleRowStatistics& sample_statistics() {
		return m_sample_statistics ;
	}
	
	bool ignore_warnings() const {
		return m_ignore_warnings ;
	}

	virtual OstreamTee& logger() {
		return m_logger ;
	}
	
	std::vector< std::size_t >& snp_filter_failure_counts() { return m_snp_filter_failure_counts ; }
	std::vector< std::size_t >& sample_filter_failure_counts() { return m_sample_filter_failure_counts ; }
	
	void write_start_banner() {
		m_logger << "\nWelcome to qc-tool\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}

	void write_end_banner() {
		m_logger << "\nThank you for using qc-tool.\n" ;
	}
		
	void write_preamble() {
		m_logger << std::setw(30) << "Input SAMPLE file:"
			<< "  \"" << m_options.input_sample_filename() << "\".\n" ;
		m_logger << std::setw(30) << "Output SAMPLE file:"
			<< "  \"" << m_options.output_sample_filename() << "\".\n" ;
		m_logger << std::setw(30) << "Sample statistic output file:"
			<< "  \"" << m_options.output_sample_stats_filename() << "\".\n" ;
		m_logger << "\n" ;

		m_logger << std::setw(30) << "Input GEN file(s):" ;
		for( std::size_t i = 0; i < m_options.gen_filename_mapper().input_files().size(); ++i ) {
			if( i > 0 ) {
				m_logger << std::string( 30, ' ' ) ;
			}
			m_logger << "  (" << std::setw(6) << m_snp_data_source->number_of_snps_in_source( i ) << " snps)  " ;
			m_logger << "\"" << m_options.gen_filename_mapper().input_files()[i] << "\"" ;				
			m_logger << "\n" ;
		}
		m_logger << std::string( 30, ' ' ) << "  (total " << m_snp_data_source->total_number_of_snps() << " snps in input).\n" ;
		if( m_options.gen_filename_mapper().input_files().size() > 1 ) {
			m_logger << "\n" ;
		}

		m_logger << std::setw(30) << "Output GEN file(s):" ;
		if( m_options.gen_filename_mapper().output_filenames().empty() ) {
			m_logger << "  (n/a)\n" ;
		}
		else {
			for( std::size_t i = 0; i < m_options.gen_filename_mapper().output_filenames().size(); ++i ) {
				if( i > 0 ) {
					m_logger << "\n" << std::string( 30, ' ' ) ;
				}
				m_logger << "  \"" << m_options.gen_filename_mapper().output_filenames()[i] << "\"" ;				
			}
			m_logger << "\n" ;
		}

		m_logger << std::setw(30) << "SNP statistic output file(s):" ;
		for( std::size_t i = 0; i < m_options.snp_stats_filename_mapper().output_filenames().size(); ++i ) {
			if( i > 0 ) {
				m_logger << std::string( 30, ' ' ) ;
			}
			m_logger << "  \"" << m_options.snp_stats_filename_mapper().output_filenames()[i] << "\"\n" ;
		}
		m_logger << "\n" ;
		m_logger << std::setw(30) << "Sample filter:" 
			<< "  " << *m_sample_filter << ".\n" ;
		m_logger << std::setw(30) << "SNP filter:"
			<< "  " << *m_snp_filter << ".\n" ;
		m_logger << "\n" ;

		m_logger << std::setw(30) << "# of samples in input files:"
			<< "  " << m_snp_data_source->number_of_samples() << ".\n" ;
		m_logger << std::setw(30) << "# of samples after filtering:"
			<< "  " << m_snp_data_source->number_of_samples() - m_indices_of_filtered_out_samples.size()
			<< " (" << m_indices_of_filtered_out_samples.size()
			<< " filtered out).\n" ;
	}
	
	void write_postamble() {
		m_logger << std::string( 72, '=' ) << "\n\n" ;
		try {
			if( m_backup_creator.backed_up_files().size() > 0 ) {
				m_logger << std::setw(36) << "I took backups of the following files before overwriting:\n" ;
				std::size_t max_length = 0u ;
				for(
					std::map< std::string, std::string >::const_iterator i = m_backup_creator.backed_up_files().begin() ;
					i != m_backup_creator.backed_up_files().end() ;
					++i
				) {
					max_length = std::max( max_length, i->first.size() ) ;
				}

				for(
					std::map< std::string, std::string >::const_iterator i = m_backup_creator.backed_up_files().begin() ;
					i != m_backup_creator.backed_up_files().end() ;
					++i
				) {
					m_logger << "  " << std::setw( max_length + 2 ) << std::left << ("\"" + i->first + "\"") << " to \"" << i->second << "\"\n" ;
				}

				m_logger << "\n" ;
				m_logger << std::string( 72, '=' ) << "\n\n" ;
			}

			m_logger << std::setw(36) << "Number of SNPs in input file(s):"
				<< "  " << m_snp_data_source->total_number_of_snps() << ".\n" ;
			if( m_snp_filter->number_of_subconditions() > 0 ) {
				for( std::size_t i = 0; i < m_snp_filter->number_of_subconditions(); ++i ) {
					m_logger << std::setw(36) << ("...which failed \"" + to_string( m_snp_filter->subcondition( i )) + "\":")
						<< "  " << m_snp_filter_failure_counts[i] << ".\n" ;
				}

				m_logger << std::setw(36) << "(total failures:"
					<< "  " << m_snp_data_source->total_number_of_snps() - m_fltrd_out_snp_data_sink->number_of_snps_written() << ").\n" ;
			}
			
			m_logger << "\n" ;

			m_logger << std::setw(36) << "Number of samples in input file(s):"
				<< "  " << m_snp_data_source->number_of_samples() << ".\n" ;
			if( m_sample_filter->number_of_subconditions() > 0 ) {
				for( std::size_t i = 0 ; i < m_sample_filter_failure_counts.size(); ++ i ) {
					m_logger << std::setw(36) << ("...which failed \"" + to_string( m_sample_filter->subcondition( i )) + "\":")
						<< "  " << m_sample_filter_failure_counts[i] << ".\n" ;
				}
				m_logger << std::setw(36) << "(total failures:" ;
				{
					std::size_t total_failures = 0 ;
					total_failures = std::accumulate( m_sample_filter_failure_counts.begin(), m_sample_filter_failure_counts.end(), total_failures ) ;
					m_logger << "  " << total_failures << ").\n" ;
				}
			}

			m_logger << "\n" ;

			if( m_options.gen_filename_mapper().output_filenames().size() > 0 ) {
				m_logger << std::setw(36) << "Output GEN files:" ;
				for( std::size_t i = 0; i < m_options.gen_filename_mapper().output_filenames().size(); ++i ) {
					if( i > 0 ) {
						m_logger << std::string( 36, ' ' ) ;
					}
					if( m_fltrd_in_snp_data_sink.get() ) {
						m_logger << "  (" << std::setw(6) << m_fltrd_in_snp_data_sink->sink(i).number_of_snps_written() << " snps)  " ;
					}
					m_logger << "\"" << m_options.gen_filename_mapper().output_filenames()[i] << "\"\n" ;
				}
				m_logger << std::string( 36, ' ' ) << "  (total " << m_fltrd_in_snp_data_sink->number_of_snps_written() << " snps).\n" ;
			}
			if( m_options.output_sample_filename() != "" ) {
				m_logger << std::setw(36) << "Output SAMPLE files:"
					<< "  \"" << m_options.output_sample_filename() << "\""
					<< "  (" << m_sample_rows.size() << " samples)\n" ;
			}
			if( m_options.snp_stats_filename_mapper().output_filenames().size() > 0 ) {
				m_logger << std::setw(36) << "SNP statistic output file(s):" ;
				for( std::size_t i = 0; i < m_options.snp_stats_filename_mapper().output_filenames().size(); ++i ) {
					if( i > 0 ) {
						m_logger << std::string( 30, ' ' ) ;
					}
					m_logger << "  \"" << m_options.snp_stats_filename_mapper().output_filenames()[i] << "\"\n" ;
				}
			}
			if( m_options.output_sample_stats_filename() != "" ) {
				m_logger << std::setw(36) << "Sample statistic output file:"
					<< "  \"" << m_options.output_sample_stats_filename() << "\".\n" ;
			}

			m_logger[ "screen" ] << std::setw( 36 ) << "\nMore details are in the log file:"
				<< "  \"" << m_options.log_filename() << "\".\n" ;
			m_logger << std::string( 72, '=' ) << "\n\n" ;
		}
		catch (...) {
			throw ;
		}
	}

	void print_progress_if_necessary() {
		double time_now = m_timer.elapsed() ;
		if( (time_now - m_last_timestamp > 1.0) || ( snp_data_source().number_of_snps_read() == snp_data_source().total_number_of_snps()) ) {
			print_progress() ;
			timestamp() ;
		}
	}

	void print_progress() {
		m_logger[ "screen" ]
			<< "\r"
			<< get_progress_bar( 30, static_cast< double >( m_snp_data_source->number_of_snps_read() ) / m_snp_data_source->total_number_of_snps() )
			<< " (" << m_snp_data_source->number_of_snps_read() << " / " << m_snp_data_source->total_number_of_snps()
			<< ", " << std::fixed << std::setprecision(1) << m_timer.elapsed() << "s)"
			<< std::flush ;
	}
	
	void timestamp() { m_last_timestamp = m_timer.elapsed() ; }

private:
	
	void setup() {
		open_log() ;
		try {
			open_gen_row_source() ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_logger << "The GEN files specified did not all have the same sample size.\n" ;
			throw ;
		}
		construct_snp_statistics() ;
		construct_sample_statistics() ;
		construct_snp_filter() ;
		construct_sample_filter() ;
		process_other_options() ;
		load_sample_rows() ;
		check_for_errors_and_warnings() ;
	}
	
	void open_log() {
		m_logger.add_stream( "screen", std::cout ) ;
		m_backup_creator.backup_file_if_necessary( m_options.log_filename() ) ;
		m_log.open( m_options.log_filename().c_str() ) ;
		if( m_log.is_open() ) {
			// Make all output go to the log as well.
			m_logger.add_stream( "log", m_log ) ;
		}
	}

	void open_gen_row_source() {
		Timer timer ;
		
		std::auto_ptr< genfile::SNPDataSourceChain > chain( new genfile::SNPDataSourceChain() ) ;
		chain->set_moved_to_next_source_callback( boost::bind( &QCToolCmdLineContext::move_to_next_output_file, this, _1 )) ;

		for( std::size_t i = 0; i < m_options.gen_filename_mapper().input_files().size(); ++i ) {
			Timer file_timer ;
			m_logger << "(Opening gen file \"" << m_options.gen_filename_mapper().input_files()[i] << "\"...)" << std::flush ;
			try {
				chain->add_source( genfile::SNPDataSource::create( m_options.gen_filename_mapper().input_files()[i] ) ) ;
			}
			catch ( genfile::FileHasTwoConsecutiveNewlinesError const& e ) {
				std::cerr << "\n!!ERROR: a GEN file was specified having two consecutive newlines.\n"
					<< "!! NOTE: popular editors, such as vim and nano, automatically add an extra newline to the file (which you can't see).\n"
					<< "!!     : Please check that each SNP in the file is terminated by a single newline.\n" ;
				throw ;
			}
			m_logger << " (" << file_timer.elapsed() << "s)\n" ;
		}

		m_snp_data_source = chain ;

		if( timer.elapsed() > 1.0 ) {
			m_logger << "Opened " << m_options.gen_filename_mapper().input_files().size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
		}
	}

	void move_to_next_output_file( std::size_t index ) {
		if( index < m_options.gen_filename_mapper().input_files().size() ) {
			if( m_options.gen_filename_mapper().output_filenames().size() > 0 ) {
				if( m_options.gen_filename_mapper().index_of_filename_corresponding_to( index ) != m_fltrd_in_snp_data_sink->index_of_current_sink() ) {
					m_fltrd_in_snp_data_sink->move_to_next_sink() ;
				}
			}
			
			if( m_options.snp_stats_filename_mapper().output_filenames().size() > 0 ) {
				if( m_options.snp_stats_filename_mapper().index_of_filename_corresponding_to( index ) != m_current_snp_stats_filename_index ) {
					open_snp_stats_file( ++m_current_snp_stats_filename_index, m_snp_statistics ) ;
				}
			}
		}
	}

	void open_snp_data_sink() {
		reset_filtered_in_snp_data_sink() ;
		open_filtered_in_snp_data_sink() ;
		reset_filtered_out_snp_data_sink() ;
		open_filtered_out_snp_data_sink() ;
	}

	void reset_filtered_in_snp_data_sink() {
		m_fltrd_in_snp_data_sink = std::auto_ptr< genfile::SNPDataSinkChain >( new genfile::SNPDataSinkChain() ) ;
	}

	void open_filtered_in_snp_data_sink() {
		if( m_options.gen_filename_mapper().output_filenames().size() == 0 ) {
			m_fltrd_in_snp_data_sink->add_sink( std::auto_ptr< genfile::SNPDataSink >( new genfile::TrivialSNPDataSink() )) ;
		}
		else {
			for( std::size_t i = 0; i < m_options.gen_filename_mapper().output_filenames().size(); ++i ) {
				std::string const& filename = m_options.gen_filename_mapper().output_filenames()[i] ;
				m_fltrd_in_snp_data_sink->add_sink( genfile::SNPDataSink::create( filename )) ;
			}
		}
	}

	void reset_filtered_out_snp_data_sink() {
		m_fltrd_out_snp_data_sink.reset( new genfile::SNPDataSinkChain() ) ;
	}

	void open_filtered_out_snp_data_sink() {
		m_fltrd_out_snp_data_sink->add_sink( std::auto_ptr< genfile::SNPDataSink >( new genfile::TrivialSNPDataSink() )) ;
	}

	void open_sample_row_source() {
		m_sample_source.reset( new NullObjectSource< SampleRow >()) ;
		if( m_options.input_sample_filename() != "" ) {
			m_sample_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_options.input_sample_filename() ))) ;
		}
	} ;

	void open_sample_row_sink() {
		m_fltrd_in_sample_sink.reset( new NullObjectSink< SampleRow >() ) ;
		if( m_options.output_sample_filename() != "" ) {
			m_fltrd_in_sample_sink.reset( new SampleOutputFile< SimpleFileObjectSink< SampleRow > >( open_file_for_output( m_options.output_sample_filename() ))) ;
		}
	}

	void open_sample_stats_sink() {
		if( m_options.output_sample_stats_filename() != "" ) {
			m_sample_stats_file = open_file_for_output( m_options.output_sample_stats_filename() ) ;
			m_sample_statistics.format_column_headers( *m_sample_stats_file ) ;
			(*m_sample_stats_file) << "\n" ;
		}
	}
	
	void reset_sample_stats_sink() {
		m_sample_stats_file.reset() ;
	}

	void open_snp_stats_file( std::size_t index, GenRowStatistics const& snp_statistics ) {
		assert( index < m_options.snp_stats_filename_mapper().output_filenames().size()) ;
		m_current_snp_stats_filename_index = index ;
		m_snp_stats_file = open_file_for_output( m_options.snp_stats_filename_mapper().output_filenames()[ index ] ) ;
		*m_snp_stats_file << "        " ;
		snp_statistics.format_column_headers( *m_snp_stats_file ) << "\n";
	}

	void reset_snp_stats_sink() {
		m_snp_stats_file.reset() ;
	}

	void construct_snp_statistics() {
		GenRowStatisticFactory::add_statistics( m_options.row_statistics_specs(), m_snp_statistics ) ;
	}

	void construct_sample_statistics() {
		std::vector< std::string > sample_statistics_specs = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-sample-stats-columns" ), "," ) ;
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
			std::vector< std::string > filenames = m_options.get_values< std::string >( "-snp-incl-list" ) ;
			std::auto_ptr< RowCondition > snp_incl_condition( new SNPInListCondition( filenames )) ;
			snp_filter->add_subcondition( snp_incl_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "-snp-excl-list" ) ) {
			std::vector< std::string > filenames = m_options.get_values< std::string >( "-snp-excl-list" ) ;
			std::auto_ptr< RowCondition > snp_incl_condition( new SNPInListCondition( filenames )) ;
			std::auto_ptr< RowCondition > snp_excl_condition( new NotRowCondition( snp_incl_condition )) ;
			snp_filter->add_subcondition( snp_excl_condition ) ;
		}
		
		m_snp_filter = snp_filter ;
		m_snp_filter_failure_counts.resize( m_snp_filter->number_of_subconditions(), 0 ) ;
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
		m_sample_filter_failure_counts.resize( m_sample_filter->number_of_subconditions(), 0 ) ;
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
	
	void load_sample_rows() {
		try {
			unsafe_load_sample_rows() ;
		}
		catch( ConditionValueNotFoundException const& ) {
			m_logger << "\n\n!! ERROR: The input sample file must contain entries for all values used to filter on.\n"
				<< "!! This includes \"missing\" and \"heterozygosity\".\n" ;
			throw ;
		}
	}
	
	void unsafe_load_sample_rows() {
		if( m_options.input_sample_filename() != "" ) {
			Timer timer ;
			m_logger << "(Preprocessing sample file...)" ;
			open_sample_row_source() ;
			SampleRow sample_row ;
			while((*m_sample_source) >> sample_row ) {
				m_sample_rows.push_back( sample_row ) ;
				if( !m_sample_filter->check_if_satisfied( sample_row )) {
					m_indices_of_filtered_out_samples.push_back( m_sample_rows.size() - 1 ) ;
				}
			}

			m_logger << "(" << std::fixed << std::setprecision(1) << timer.elapsed() << "s)\n" ;
			if( m_sample_rows.size() != m_snp_data_source->number_of_samples() ) {
				throw NumberOfSamplesMismatchException() ;
			}
		}
		else {
			m_sample_rows.resize( m_snp_data_source->number_of_samples() ) ;
		}
	}
	
	void check_for_errors_and_warnings() {
		check_for_errors() ;
		check_for_warnings() ;
	}

	void check_for_errors() {
		if( m_options.gen_filename_mapper().input_files().size() == 0 ) {
			m_errors.push_back( "At least one GEN input file must be supplied." ) ;
		}

		for( std::size_t i = 0; i < m_options.gen_filename_mapper().input_files().size(); ++i ) {
			for( std::size_t j = 0; j < m_options.gen_filename_mapper().output_filenames().size(); ++j ) {
				if( strings_are_nonempty_and_equal( m_options.gen_filename_mapper().output_filenames()[j], m_options.gen_filename_mapper().input_files()[i] )) {
					m_errors.push_back( "Output GEN file \"" + m_options.gen_filename_mapper().output_filenames()[j] +"\" also specified as input GEN file." ) ;
					break ;
				}
			}
			if( strings_are_nonempty_and_equal( m_options.output_sample_filename(), m_options.gen_filename_mapper().input_files()[i] )) {
				m_errors.push_back( "Output SAMPLE file \"" + m_options.output_sample_filename() +"\" also specified as input GEN file." ) ;
				break ;
			}
			for( std::size_t j = 0; j < m_options.snp_stats_filename_mapper().output_filenames().size(); ++j ) {
				if( strings_are_nonempty_and_equal( m_options.snp_stats_filename_mapper().output_filenames()[j], m_options.gen_filename_mapper().input_files()[i] )) {
					m_errors.push_back( "Output GEN statistic file \"" + m_options.snp_stats_filename_mapper().output_filenames()[j] +"\" also specified as input GEN file." ) ;
					break ;
				}
			}
			if( strings_are_nonempty_and_equal( m_options.output_sample_stats_filename(), m_options.gen_filename_mapper().input_files()[i] )) {
				m_errors.push_back( "Output SAMPLE statistic file \"" + m_options.output_sample_stats_filename() +"\" also specified as input GEN file." ) ;
				break ;
			}
		}
		for( std::size_t j = 0; j < m_options.gen_filename_mapper().output_filenames().size(); ++j ) {
			if( strings_are_nonempty_and_equal( m_options.output_sample_filename(), m_options.gen_filename_mapper().output_filenames()[j] )) {
				m_errors.push_back( "The GEN and SAMPLE output filenames must differ." ) ;
			}
		}
		for( std::size_t j = 0; j < m_options.snp_stats_filename_mapper().output_filenames().size(); ++j ) {
			if( strings_are_nonempty_and_equal( m_options.snp_stats_filename_mapper().output_filenames()[j], m_options.output_sample_stats_filename() )) {
				m_errors.push_back( "The gen statistic and sample statistic filenames must differ." ) ;
			}
		}
		if( m_options.input_sample_filename() == "" && m_sample_filter->number_of_subconditions() != 0 ) {
			m_errors.push_back( "To filter on samples, please supply an input sample file." ) ;
		}
	}
	
	void check_for_warnings() {
		if( (m_options.output_sample_filename() != "" || m_options.output_sample_stats_filename() != "") && m_options.gen_filename_mapper().input_files().size() != 23 ) {
			m_warnings.push_back( "You are outputting a sample or sample statistic file, but the number of gen files is not 23.\n"
			"   (I suspect there is not the whole genomes' worth of data?)" ) ;
		}
		if( m_options.output_sample_stats_filename() != "" && m_options.input_sample_filename() == "" ) {
			m_warnings.push_back( "You are outputting a sample statistic file, but no input sample file has been supplied.\n"
			"   Statistics will be output but the ID fields will be left blank.") ;
		}
		if( m_options.gen_filename_mapper().output_filenames().size() == 0 && m_options.output_sample_filename() == "" && m_options.snp_stats_filename_mapper().output_filenames().size() == 0 && m_options.output_sample_stats_filename() == "" ) {
			m_warnings.push_back( "You have not specified any output files.  This will produce only console output." ) ;
		}
		if( (m_options.gen_filename_mapper().output_filenames().size() > 0) &&  (m_snp_filter->number_of_subconditions() == 0) && (m_sample_filter->number_of_subconditions() == 0)) {
			m_warnings.push_back( "You have specified output GEN files, but no filters.\n"
			 	"The output GEN files will contain the same data as the input ones.") ;
		}
		if(m_options.gen_filename_mapper().output_filenames().size() == 0 && ((m_snp_filter->number_of_subconditions() > 0) || (m_sample_filter->number_of_subconditions() > 0))) {
			m_warnings.push_back( "You have not specified filters, but no output GEN files." ) ;
		}
		if(m_options.output_sample_filename() == "" && (m_sample_filter->number_of_subconditions() > 0)) {
			m_warnings.push_back( "You have not specified sample filters, but no output sample file.\n" ) ;
		}
		if( strings_are_nonempty_and_equal( m_options.output_sample_filename(), m_options.input_sample_filename() )) {
			m_warnings.push_back( "The input sample file \"" + m_options.input_sample_filename() + "\" will be overwritten with the output sample file.\n"
				"  (A backup will be taken, but if overwriting is not desired, please use the -os option to choose a different filename.)" ) ;
		}
	}
	
	bool strings_are_nonempty_and_equal( std::string const& left, std::string const& right ) {
		return (!left.empty()) && (!right.empty()) && (left == right) ;
	}

private:
	Timer m_timer ;
	double m_last_timestamp ;
	
	QCToolOptionProcessor m_options ;

	std::ofstream m_log ;
	OstreamTee m_logger ;
	
	std::auto_ptr< genfile::SNPDataSourceChain > m_snp_data_source ;
	std::auto_ptr< genfile::SNPDataSinkChain > m_fltrd_in_snp_data_sink ;
	std::auto_ptr< genfile::SNPDataSinkChain > m_fltrd_out_snp_data_sink ;
	std::auto_ptr< ObjectSource< SampleRow > > m_sample_source ;

	std::vector< SampleRow > m_sample_rows ;

	std::auto_ptr< ObjectSink< SampleRow > > m_fltrd_in_sample_sink ;
	std::auto_ptr< ObjectSink< SampleRow > > m_fltrd_out_sample_sink ;
	OUTPUT_FILE_PTR m_snp_stats_file ;
	OUTPUT_FILE_PTR m_sample_stats_file ;

	std::vector< std::size_t > m_snp_filter_failure_counts ;
	std::vector< std::size_t > m_sample_filter_failure_counts ;
	
	bool m_ignore_warnings ; 
	
	std::size_t m_current_snp_stats_filename_index ;
	
	std::auto_ptr< AndRowCondition > m_snp_filter ;
	std::auto_ptr< AndRowCondition > m_sample_filter ;

	GenRowStatistics m_snp_statistics ;
	SampleRowStatistics m_sample_statistics ;
	
	std::vector< std::size_t > m_indices_of_filtered_out_samples ;
	
	FileBackupCreator m_backup_creator ;
	
	std::vector< std::string > m_warnings ;
	std::vector< std::string > m_errors ;
} ;


struct QCToolProcessor
{
public:
	QCToolProcessor( QCToolContext & context )
		: m_context( context )
	{
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
		process_gen_rows() ;
		process_sample_rows() ;
		construct_plots() ;
	}

	void process_gen_rows() {
		m_context.logger() << "Processing SNPs...\n" ;
		Timer timer ;

		InternalStorageGenRow row ;
		std::size_t number_of_snps_processed = 0 ;
		for(
			row.set_number_of_samples( m_context.snp_data_source().number_of_samples() ) ;
			row.read_from_source( m_context.snp_data_source() ) ;
			row.set_number_of_samples( m_context.snp_data_source().number_of_samples() )
		) {
			preprocess_gen_row( row ) ;
			process_gen_row( row, ++number_of_snps_processed ) ;
			accumulate_per_column_amounts( row, m_per_column_amounts ) ;
			m_context.print_progress_if_necessary() ;
		}
		m_context.logger() << "\n" ;
		assert( m_context.snp_data_source().number_of_snps_read() == m_context.snp_data_source().total_number_of_snps() ) ;
		m_context.logger() << "Processed " << m_context.snp_data_source().total_number_of_snps() << " SNPs in "
			<< std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
		if( m_context.snp_filter().number_of_subconditions() > 0 ) {
			m_context.logger() << "(" << m_context.fltrd_in_snp_data_sink().number_of_snps_written() << " of " << m_context.snp_data_source().total_number_of_snps() << " SNPs passed the filter.)\n" ;
		}
	}
	
	void preprocess_gen_row( InternalStorageGenRow& row ) const {
		check_gen_row( row ) ;
		row.filter_out_samples_with_indices( m_indices_of_filtered_out_samples ) ;
	}
	
	void check_gen_row( GenRow& row ) const {
		check_gen_row_has_correct_number_of_samples( row ) ;
	}

	void check_gen_row_has_correct_number_of_samples( GenRow& row ) const {
		if( row.number_of_samples() != m_context.sample_rows().size() ) {
			throw GenAndSampleFileMismatchException() ;
		}
	}

	void process_gen_row( GenRow const& row, std::size_t row_number ) {
		m_context.snp_statistics().process( row ) ;
		if( m_context.snp_filter().check_if_satisfied( m_context.snp_statistics() )) {
			row.write_to_sink( m_context.fltrd_in_snp_data_sink() ) ;
			output_gen_row_stats( m_context.snp_statistics(), row_number ) ;
		}
		else {
			row.write_to_sink( m_context.fltrd_out_snp_data_sink() ) ;
			do_snp_filter_diagnostics( m_context.snp_statistics(), row_number ) ;
		}
	}

	void output_gen_row_stats( GenotypeAssayStatistics const& row_statistics, std::size_t row_number ) {
		if( m_context.snp_statistics().size() > 0 ) {
			m_context.snp_stats_sink()
				<< std::setw(8) << std::left << row_number
				<< m_context.snp_statistics() << "\n" ;
		}
	}

	void do_snp_filter_diagnostics( GenRowStatistics const& row_statistics, std::size_t const row_index ) {
		std::ostream& log = m_context.logger()["log"] ;
		log
			<< "Filtered out snp #" << row_index << " (" << row_statistics.row().SNPID() << " " << row_statistics.row().RSID() << " " << row_statistics.row().SNP_position() << ")"
			<< " because it does not satisfy " ;
		for( std::size_t i = 0, failed_condition_count = 0; i < m_context.snp_filter().number_of_subconditions() ; ++i ) {
			if( !m_context.snp_filter().subcondition(i).check_if_satisfied( row_statistics )) {

				++m_context.snp_filter_failure_counts()[ i ] ;

				if( failed_condition_count > 0 ) {
					log << " or " ;
				}
				log << "\"" << m_context.snp_filter().subcondition(i) << "\"" ;
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
	

	void process_sample_rows() {
		m_context.logger() << "Processing samples...\n" ;
		Timer timer ;

		apply_sample_filter() ;
		assert( m_context.sample_rows().size() == m_per_column_amounts.size() ) ;

		for( std::size_t i = 0 ; i < m_per_column_amounts.size(); ++i ) {
			SampleRow& sample_row = m_context.sample_rows()[i] ;
			m_context.sample_statistics().process( sample_row, m_per_column_amounts[i], m_context.snp_data_source().total_number_of_snps() ) ;
			m_context.sample_stats_sink() << std::setw(8) << (i+1) << m_context.sample_statistics() << "\n" ;
			m_context.sample_statistics().add_to_sample_row( sample_row, "missing" ) ;
			m_context.sample_statistics().add_to_sample_row( sample_row, "heterozygosity" ) ;
			m_context.fltrd_in_sample_sink() << sample_row ;
		}

		m_context.logger() << "Processed " << m_context.snp_data_source().number_of_samples() << " samples in " << std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
		if( m_context.sample_filter().number_of_subconditions() > 0 ) {
			m_context.logger() << "(" << m_context.sample_rows().size() << " of " << m_context.snp_data_source().number_of_samples() << " samples passed the filter.)\n" ;
		}
		m_context.logger() << "\n" ;
	}

	void apply_sample_filter() {
		for( std::size_t pre_filter_i = 0, post_filter_i = 0; post_filter_i < m_context.sample_rows().size(); ++pre_filter_i ) {
			if( sample_row_is_filtered_out( pre_filter_i ) ) {
				do_sample_filter_diagnostics( m_context.sample_rows()[ post_filter_i ], pre_filter_i ) ;
				m_context.sample_rows().erase( m_context.sample_rows().begin() + post_filter_i ) ;
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
		std::ostream& log = m_context.logger()["log"] ;
		log
			<< "Filtered out sample row " << sample_row_index << " (" << sample_row.ID1() << " " << sample_row.ID2() << ")"
			<< " because it does not satisfy " ;
		for( std::size_t i = 0, failed_condition_count = 0; i < m_context.sample_filter().number_of_subconditions() ; ++i ) {
			if( !m_context.sample_filter().subcondition(i).check_if_satisfied( sample_row )) {
				++m_context.sample_filter_failure_counts()[ i ] ;
				if( failed_condition_count > 0 ) {
					log << " or " ;
				}
				log << "\"" << m_context.sample_filter().subcondition(i) << "\"" ;
				++failed_condition_count ;
			}
		}
		log << ".\n" ;
	}

	void construct_plots() {
		// not implemented.
	}

private:
	QCToolContext& m_context ;

	std::size_t m_number_of_sample_file_rows ;
	std::size_t m_number_of_filtered_in_snps ;
	std::vector< GenotypeProportions > m_per_column_amounts ;
	std::vector< std::size_t > m_indices_of_filtered_out_samples ;
} ;


int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		QCToolCmdLineContext context( argc, argv ) ;
		QCToolProcessor processor( context ) ;
		processor.process() ;
    }
	catch( HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
