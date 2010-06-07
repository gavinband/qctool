/*
 * This program, qctool, selects rows from a GEN file according to certain criteria.
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
#include "SNPIDMatchesCondition.hpp"
#include "SampleInListCondition.hpp"
#include "FileUtil.hpp"
#include "GenRowStatistics.hpp"
#include "SampleRowStatistics.hpp"
#include "ObjectSource.hpp"
#include "SimpleFileObjectSource.hpp"
#include "SimpleFileObjectSink.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/SNPDataSinkChain.hpp"
#include "genfile/TrivialSNPDataSink.hpp"
#include "genfile/CategoricalCohortIndividualSource.hpp"

#include "SampleOutputFile.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "wildcard.hpp"
#include "string_utils/string_utils.hpp"
#include "string_utils/parse_utils.hpp"
#include "progress_bar.hpp"
#include "FileBackupCreator.hpp"
#include "InputToOutputFilenameMapper.hpp"
#include "OstreamTee.hpp"
#include "null_ostream.hpp"
#include "SNPPositionSink.hpp"
#include "SNPIDSink.hpp"
#include "statfile/RFormatStatSink.hpp"
#include "SampleIDSink.hpp"
#include "QCToolContext.hpp"
#include "QCToolProcessor.hpp"


namespace globals {
	std::string const program_name = "qctool" ;
}

struct NumberOfSamplesMismatchException: public QCToolProcessorException
{
	char const* what() const throw() {return "NumberOfSamplesMismatchException" ; }
} ;

// Thrown to indicate that the numbers of input and output files specified on the command line differ.
struct QCToolFileCountMismatchError: public QCToolProcessorException
{
	char const* what() const throw() {return "QCToolFileCountMismatchError" ; }
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
								"matches numbers from 1 to 100.  For example, \"qctool -g myfile_#.gen\" will find all files of "
								"the form \"myfile_N.gen\", where N can be 1, 002, 099, etc." )
	        .set_is_required()
			.set_takes_values()
			.set_maximum_number_of_repeats(23) ;

	    options[ "-og" ]
	        .set_description( 	"Override the auto-generated path(s) of the output gen file(s) to use when filtering.  "
								"(By default, the paths are formed by adding \".fltrd\" to the input gen filename(s).)  "
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
			.set_description( "Calculate and output snp-wise statistics." ) ;
	    options[ "-snp-stats-file" ]
	        .set_description( 	"Override the auto-generated path(s) of the snp-stats file to use when outputting snp-wise statistics.  "
								"(By default, the paths are formed by adding \".snp-stats\" to the input gen filename(s).)  "
								"If used, this option must appear as many times as the -g option.  "
	 							"If the corresponding occurence of -g uses a '#' wildcard character, the '#' character can "
								"also be used here to specify numbered output files corresponding to the input files." )
	        .set_takes_values()
			.set_maximum_number_of_repeats(23) ;

		options[ "-snp-stats-columns" ]
	        .set_description( "Comma-seperated list of columns to output in the snp-wise statistics file.  "
	 						"By default, the columns are: "
							"SNPID, RSID, position, minor_allele, major_allele, MAF, HWE, missing, and information" )
			.set_takes_single_value()
			.set_default_value( "SNPID, RSID, position, minor_allele, major_allele, AA, AB, BB, MAF, HWE, missing, information" ) ;

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
			.set_description( "Filter out SNPs with -log10( HWE p-value ) greater than or equal to the value specified.")
			.set_takes_single_value() ;
		options[ "-info" ]
			.set_description( "Filter out SNPs with Fisher information lying outside the given range.")
			.set_number_of_values_per_use( 2 ) ;
		options[ "-snp-missing-rate" ]
			.set_description( "Filter out SNPs with missing data rate greater than or equal to the value specified.")
			.set_takes_single_value() ;
		options[ "-snp-interval" ]
			.set_description( "Filter out SNPs with position outside the interval [a,b]." )
			.set_number_of_values_per_use( 2 ) ;
		options[ "-maf" ]
			.set_description( "Filter out SNPs whose minor allele frequency lies outside the interval [a,b]." )
			.set_number_of_values_per_use( 2 ) ;
		options[ "-snpid-filter" ]
			.set_description( "Filter out snps whose SNPID doesn't match the given argument.  "
						"The argument must be a string, optionally containing a wildcard character ('*')"
						" which matches any sequence of characters.")
			.set_takes_single_value() ;

		// Sample filtering options
		options.declare_group( "Sample filtering options" ) ;
		options[ "-sample-missing-rate" ]
			.set_description( "Filter out samples with missing data rate greater than the value specified.")
			.set_takes_single_value() ;
		options[ "-heterozygosity" ]
			.set_description( "Filter out samples with heterozygosity outside the inteval [a,b]." )
			.set_number_of_values_per_use( 2 ) ;

		// Inclusion / exclusion list options
		options.declare_group( "Inclusion / exclusion list options" ) ;
		options[ "-sample-incl-list"]
			.set_description( "Filter out samples whose sample ID does not lie in the given file.")
			.set_takes_single_value() ;
		options[ "-sample-excl-list"]
			.set_description( "Filter out samples whose sample ID lies in the given file.")
			.set_takes_single_value() ;
		options[ "-write-sample-excl-list" ]
			.set_description( "Do not apply sample filters directly.  Instead, write a file containing a list of the ids"
			"  of individuals which would be filtered out by the filter." ) ;
		options[ "-write-sample-excl-list-file" ]
			.set_description( "Override default name of the file to use in -write-sample-excl-list" )
			.set_takes_single_value() ;
		options[ "-snp-incl-list" ]
			.set_description( "Filter out SNPs whose SNP ID or RSID does not lie in the given file(s).")
			.set_takes_values() 
			.set_number_of_values_per_use( 1 )
			.set_maximum_number_of_repeats( 100 ) ;
		options[ "-snp-excl-list" ]
			.set_description( "Filter out SNPs whose SNP ID or RSID lies in the given file(s).")
			.set_takes_values() 
			.set_number_of_values_per_use( 1 )
			.set_maximum_number_of_repeats( 100 ) ;
		options [ "-write-snp-excl-list" ]
			.set_description( "Don't apply the filter; instead, write files containing the positions of SNPs that would be filtered out."
			"  These files are suitable for use as the input to -snp-incl-list option on a subsequent run." ) ;
		options [ "-write-snp-excl-list-file" ]
	        .set_description( 	"Override the auto-generated path(s) of the file to use when outputting the positions of filtered out SNPs.  "
								"(By default, the paths are formed by adding \".snp-excl-list\" to the input gen filename(s).)  "
								"If used, this option must appear as many times as the -g option.  "
	 							"If the corresponding occurence of -g uses a '#' wildcard character, the '#' character can "
								"also be used here to specify numbered output files corresponding to the input files." )
			.set_takes_single_value() ;
		options[ "-apply-excl-lists" ]
			.set_description( "Apply inclusion / exclusion lists supplied as arguments to the options "
			"-sample-excl-list, -sample-incl-list, -snp-excl-list, -snp-incl-list, to produce new sample / GEN files." ) ;

		// Other options
		options.declare_group( "Other options" ) ;
		options [ "-force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
		options [ "-log" ]
			.set_description( "Override the default path of the log file written by " + globals::program_name + "."
				"  By default, this is " + globals::program_name + ".log." )
			.set_takes_single_value()
			.set_default_value( globals::program_name + ".log" ) ;
		options [ "-plot" ]
			.set_description( "Path of file to produce plots in.")
			.set_takes_single_value() ;

		options.option_excludes_option( "-snp-stats", "-og" ) ;
		options.option_excludes_option( "-snp-stats", "-os" ) ;
		options.option_excludes_group( "-snp-stats", "SNP filtering options" ) ;
		options.option_excludes_group( "-sample-stats", "Sample filtering options" ) ;
		options.option_excludes_option( "-sample-stats", "-og" ) ;

		options.option_implies_option( "-write-sample-excl-list-file", "-write-sample-excl-list" ) ;
		options.option_excludes_option( "-write-sample-excl-list", "-os" ) ;
		options.option_excludes_option( "-write-sample-excl-list", "-sample-stats" ) ;

		options.option_implies_option( "-write-snp-excl-list-file", "-write-snp-excl-list" ) ;
		options.option_excludes_option( "-write-snp-excl-list", "-og" ) ;
		options.option_excludes_option( "-write-snp-excl-list", "-snp-stats" ) ;

		options.option_implies_option( "-sample-excl-list", "-s" ) ;
		options.option_implies_option( "-sample-incl-list", "-s" ) ;
		options.option_implies_option( "-snp-excl-list", "-g" ) ;
		options.option_implies_option( "-snp-incl-list", "-g" ) ;
	}
	
	void process( int argc, char** argv ) {
		CmdLineOptionProcessor::process( argc, argv ) ;
		process_filenames() ;
	}

	std::string const& input_sample_filename() const { return m_input_sample_filename ; }
	std::string const& output_sample_filename() const { return m_output_sample_filename ; }
	std::string const& output_sample_excl_list_filename() const { return m_output_sample_excl_list_filename ; }
	std::string const& output_sample_stats_filename() const { return m_sample_statistic_filename ; }
	std::string const& plot_filename() const { return m_plot_filename ; }
	std::string const log_filename() const { return get_value< std::string > ( "-log" ) ; }
	InputToOutputFilenameMapper const& gen_filename_mapper() const { return m_gen_file_mapper ; }
	InputToOutputFilenameMapper const& snp_stats_filename_mapper() const { return m_snp_stats_sink_mapper ; }
	InputToOutputFilenameMapper const& snp_excl_list_filename_mapper() const { return m_output_snp_excl_file_mapper ; }
	std::vector< std::string > row_statistics_specs() const {
		return string_utils::split_and_strip_discarding_empty_entries( get_value< std::string >( "-snp-stats-columns" ), "," ) ;
	}
	std::vector< std::string > sample_statistics_specs() const {
		return string_utils::split_and_strip_discarding_empty_entries( get_value< std::string >( "-sample-stats-columns" ), "," ) ;
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
			m_snp_stats_sink_mapper.add_filename_pair( input_gen_filenames_supplied[i], output_snp_stats_filenames_supplied[i] ) ;
			m_output_snp_excl_file_mapper.add_filename_pair( input_gen_filenames_supplied[i], output_snp_excl_filenames_supplied[i] ) ;
		}
	}

	void get_sample_related_filenames() {
		if( check_if_option_was_supplied( "-s" ) ) {
			m_input_sample_filename = get_value< std::string >( "-s" ) ;
		}

		// We need to write a sample stats file if -sample-stats was given.
		if( check_if_option_was_supplied( "-sample-stats" ) ) {
			if( check_if_option_was_supplied( "-sample-stats-file" )) {
				m_sample_statistic_filename = get_value< std::string >( "-sample-stats-file" ) ;
			}
			else {
				m_sample_statistic_filename = strip_sample_file_extension_if_present( m_input_sample_filename ) + ".sample-stats";
			}
		}
		// Otherwise, we need to write a sample exclusion list file if -write-sample-excl-list was given.
		else if( check_if_option_was_supplied( "-write-sample-excl-list" )) {
			if( check_if_option_was_supplied( "-write-sample-excl-list-file" )) {
				m_output_sample_excl_list_filename = get_value< std::string > ( "-write-sample-excl-list-file" ) ;
			}
			else {
				m_output_sample_excl_list_filename = strip_sample_file_extension_if_present( m_input_sample_filename ) + ".sample-excl-list" ;
			}
		}

		// We need to write a sample file if either
		// -write-sample-excl-list is NOT given
		// AND EITHER
		//	 * -sample-stats is given, (in which case we modify the input sample file unless overridden)
		//   * OR some sample filters are given
		//   * OR -sample-excl-list and -apply-excl-lists are both given.
		bool sample_filtering
			= check_if_option_was_supplied_in_group( "Sample filtering options" )
			|| ( check_if_option_was_supplied( "-sample-excl-list" ) && check_if_option_was_supplied( "-apply-excl-lists" )) ;

		if( !check_if_option_was_supplied( "-write-sample-excl-list" )
			&& ( check_if_option_was_supplied( "-sample-stats" ) || sample_filtering )
		) {
			if( check_if_option_was_supplied( "-os" )) {
				m_output_sample_filename = get_value< std::string >( "-os" ) ;
			}
			else {
				if( sample_filtering ) {
					m_output_sample_filename = strip_sample_file_extension_if_present( m_input_sample_filename ) + ".fltrd.sample";
				}
				else {
					m_output_sample_filename = m_input_sample_filename ;
				}
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
		// We need to produce output gen files if
		//  * -write-snp-excl-list is not set
		//  AND EITHER
		//    * some sample filters are given (but not -write-sample-excl-list)
		//    * OR some SNP filters are given
		//    * OR -sample-excl-list or -snp-excl-list ios given and the -apply-excl-lists option is given.
		if( !check_if_option_was_supplied( "-write-snp-excl-list" )
			&& (
				(check_if_option_was_supplied_in_group( "Sample filtering options" ) && !check_if_option_was_supplied( "-write-sample-excl-list" ))
				||
				check_if_option_was_supplied_in_group( "SNP filtering options" )
				||
				((check_if_option_was_supplied( "-sample-excl-list" ) || check_if_option_was_supplied( "-sample-incl-list" ) || check_if_option_was_supplied( "-snp-excl-list" ) || check_if_option_was_supplied( "-snp-incl-list" )) && check_if_option_was_supplied( "-apply-excl-lists" ))
			)
		) {
			if( check_if_option_was_supplied( "-og" ) ) {
				result = get_values< std::string >( "-og" ) ;
			} else {
				for( std::size_t i = 0; i < input_gen_filenames_supplied.size(); ++i ) {
					result[i]
						= genfile::strip_gen_file_extension_if_present( input_gen_filenames_supplied[i] )
						+ ".fltrd"
						+ genfile::get_gen_file_extension_if_present( input_gen_filenames_supplied[i] ) ;
				}
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
		if( check_if_option_was_supplied( "-write-snp-excl-list" ) ) {
			if( check_if_option_was_supplied( "-write-snp-excl-list-file" )) {
				result = get_values< std::string >( "-write-snp-excl-list-file" ) ;
			}
			else {
				for( std::size_t i = 0; i < input_gen_filenames_supplied.size() ; ++i ) {
					result[i] = genfile::strip_gen_file_extension_if_present( input_gen_filenames_supplied[i] ) + ".snp-excl-list" ;
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
	std::string m_output_sample_excl_list_filename ; 
	std::string m_sample_statistic_filename ;
	std::string m_plot_filename ;
	std::string m_log_filename ;
	InputToOutputFilenameMapper m_gen_file_mapper ;
	InputToOutputFilenameMapper m_snp_stats_sink_mapper ;
	InputToOutputFilenameMapper m_output_snp_excl_file_mapper ;
} ;


struct QCToolCmdLineContext: public QCToolContext
{
	QCToolCmdLineContext( int argc, char** argv ) {
		try {
			timestamp() ;
			m_options.process( argc, argv ) ;
			open_log() ;
			write_start_banner() ;
			setup() ;
			timestamp() ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_logger << "\nError: The GEN files specified did not all have the same sample size.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		} 
		catch ( FileError const& e ) {
			m_logger << "\nFile handling exception: " << e.what() << ": relating to file \"" << e.filename() << "\".\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch ( OptionValueInvalidException const& e ) {
			m_logger << "\nError: " << e.what() << "."
			 	<< "  (Note: " << e.option() << " takes " << e.values().size() << " values.)\n";
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch ( OptionProcessingException const& e ) {
			m_logger << "\nError: " << e.what() << ": relating to option \"" << e.option() << "\".\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
	}
	
	~QCToolCmdLineContext() {
		write_postamble() ;
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

	statfile::BuiltInTypeStatSink& snp_stats_sink() const {
		assert( m_snp_stats_sink.get() ) ;
		return *m_snp_stats_sink ;
	}

	statfile::BuiltInTypeStatSink& sample_stats_sink() const {
		assert( m_sample_stats_sink.get() ) ;
		return *m_sample_stats_sink ;
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
	
	std::vector< std::size_t > const& indices_of_filtered_out_samples() const { return m_indices_of_filtered_out_samples ; }
	
	bool ignore_warnings() const {
		return m_ignore_warnings ;
	}

	virtual OstreamTee& logger() {
		return m_logger ;
	}
	
	std::vector< std::size_t >& snp_filter_failure_counts() { return m_snp_filter_failure_counts ; }
	std::vector< std::size_t >& sample_filter_failure_counts() { return m_sample_filter_failure_counts ; }

	void write_start_banner() {
		m_logger << "\nWelcome to qctool\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}

	void write_end_banner() {
		m_logger << "\nThank you for using qctool.\n" ;
	}

	void write_preamble() {
		m_logger << std::string( 72, '=' ) << "\n\n" ;

		m_logger << std::setw(30) << "Input SAMPLE file:"
			<< "  \"" << format_filename( m_options.input_sample_filename()) << "\".\n" ;
		m_logger << std::setw(30) << "Output SAMPLE file:"
			<< "  \"" << format_filename( m_options.output_sample_filename()) << "\".\n" ;
		m_logger << std::setw(30) << "Sample statistic output file:"
			<< "  \"" << format_filename( m_options.output_sample_stats_filename()) << "\".\n" ;
		m_logger << std::setw(30) << "Sample exclusion output file:"
			<< "  \"" << format_filename( m_options.output_sample_excl_list_filename()) << "\".\n" ;
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

		m_logger << std::setw(30) << "Output SNP position file(s):" ;
		if( m_options.snp_excl_list_filename_mapper().output_filenames().empty() ) {
			m_logger << "  (n/a)\n" ;
		}
		else {
			for( std::size_t i = 0; i < m_options.snp_excl_list_filename_mapper().output_filenames().size(); ++i ) {
				if( i > 0 ) {
					m_logger << "\n" << std::string( 30, ' ' ) ;
				}
				m_logger << "  \"" << m_options.snp_excl_list_filename_mapper().output_filenames()[i] << "\"" ;				
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

		m_logger << "\n" << std::string( 72, '=' ) << "\n\n" ;

		if( !m_errors.empty() ) {
			for( std::size_t i = 0; i < m_errors.size(); ++i ) {
				m_logger << "!! ERROR: " << m_errors[i] << "\n\n" ;
			}
			m_logger << "!! Please correct the above errors and re-run qctool.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}

		if( !m_warnings.empty() ) {
			for( std::size_t i = 0; i < m_warnings.size(); ++i ) {
				m_logger << "!! WARNING: " << m_warnings[i] << "\n\n" ;
			}
			if( m_ignore_warnings ) {
				m_logger << "!! Warnings were encountered, but proceeding anyway as -force was supplied.\n" ;
				m_logger << "\n" << std::string( 72, '=' ) << "\n\n" ;
			}
			else {
				m_logger << "!! Warnings were encountered.  To proceed anyway, please run again with the -force option.\n" ;
				throw HaltProgramWithReturnCode( -1 ) ;
			}
		}
		
	}
	
	std::string format_filename( std::string const& filename ) {
		if( filename == "" ) {
			return "(n/a)" ;
		}
		else {
			return filename ;
		}
	}
	
	void write_postamble() {
		m_logger << std::string( 72, '=' ) << "\n\n" ;
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
				m_logger << std::setw(36) << ("...which failed \"" + string_utils::to_string( m_snp_filter->subcondition( i )) + "\":")
					<< "  " << m_snp_filter_failure_counts[i] << ".\n" ;
			}

			m_logger << std::setw(36) << "(total failures:"
				<< "  " << m_fltrd_out_snp_data_sink->number_of_snps_written() << ").\n" ;
		}
		
		m_logger << "\n" ;

		m_logger << std::setw(36) << "Number of samples in input file(s):"
			<< "  " << m_snp_data_source->number_of_samples() << ".\n" ;
		if( m_sample_filter->number_of_subconditions() > 0 ) {
			for( std::size_t i = 0 ; i < m_sample_filter_failure_counts.size(); ++ i ) {
				m_logger << std::setw(36) << ("...which failed \"" + string_utils::to_string( m_sample_filter->subcondition( i )) + "\":")
					<< "  " << m_sample_filter_failure_counts[i] << ".\n" ;
			}
			m_logger << std::setw(36) << "(total failures:" << "  " << m_indices_of_filtered_out_samples.size() << ").\n" ;
		}

		m_logger << "\n" ;

		if( m_options.output_sample_excl_list_filename() != "" ) {
			m_logger << std::setw(36) << "Output sample exclusion list:"
				<< " " << m_options.output_sample_excl_list_filename() << " ("
				<< m_indices_of_filtered_out_samples.size() << " samples).\n" ;
		}

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
		
		if( m_options.snp_excl_list_filename_mapper().output_filenames().size() > 0 ) {
			m_logger << std::setw(36) << "Output SNP exclusion list(s):" ;
			for( std::size_t i = 0; i < m_options.snp_excl_list_filename_mapper().output_filenames().size(); ++i ) {
				if( i > 0 ) {
					m_logger << std::string( 36, ' ' ) ;
				}
				if( m_fltrd_out_snp_data_sink.get() ) {
					m_logger << "  (" << std::setw(6) << m_fltrd_out_snp_data_sink->sink(i).number_of_snps_written() << " snps)  " ;
				}
				m_logger << "  \"" << m_options.snp_excl_list_filename_mapper().output_filenames()[i] << "\"\n" ;				
			}
			m_logger << std::string( 36, ' ' ) << "  (total " << m_fltrd_out_snp_data_sink->number_of_snps_written() << " snps).\n" ;
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
					m_logger << std::string( 36, ' ' ) ;
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
			construct_snp_statistics() ;
			construct_sample_statistics() ;
			construct_snp_filter() ;
			construct_sample_filter() ;
			process_other_options() ;

			open_snp_data_source() ;

			load_sample_rows( m_snp_data_source->number_of_samples() ) ;

			check_for_errors_and_warnings() ;

			write_preamble() ;
			open_sample_row_sink() ;
			open_snp_data_sinks() ;
			open_snp_stats_sink( 0, m_snp_statistics ) ;
			open_sample_stats_sink() ;
	}
	
	void open_log() {
		m_logger.add_stream( "screen", std::cout ) ;
		m_backup_creator.backup_file_if_necessary( m_options.log_filename() ) ;
		m_log.open( m_options.log_filename().c_str() ) ;
		if( !m_log.is_open() ) {
			std::cout << m_options.log_filename() << ".\n" ;
			throw FileNotOpenedError( m_options.log_filename() ) ;
		}
		// Make all output go to the log as well.
		m_logger.add_stream( "log", m_log ) ;
	}

	void open_snp_data_source() {
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
				if( m_options.gen_filename_mapper().filename_corresponding_to( index ) != m_fltrd_in_snp_data_sink->index_of_current_sink() ) {
					m_fltrd_in_snp_data_sink->move_to_next_sink() ;
				}
			}
			
			if( m_options.snp_stats_filename_mapper().output_filenames().size() > 0 ) {
				if( m_options.snp_stats_filename_mapper().filename_corresponding_to( index ) != m_current_snp_stats_filename_index ) {
					open_snp_stats_sink( ++m_current_snp_stats_filename_index, m_snp_statistics ) ;
				}
			}
			
			if( m_options.snp_excl_list_filename_mapper().output_filenames().size() > 0 ) {
				if( m_options.snp_excl_list_filename_mapper().filename_corresponding_to( index ) != m_fltrd_out_snp_data_sink->index_of_current_sink() ) {
					m_fltrd_out_snp_data_sink->move_to_next_sink() ;
				}
			}
		}
	}

	void open_snp_data_sinks() {
		open_filtered_in_snp_data_sink() ;
		open_filtered_out_snp_data_sink() ;
	}

	void reset_filtered_in_snp_data_sink() {
		m_fltrd_in_snp_data_sink = std::auto_ptr< genfile::SNPDataSinkChain >( new genfile::SNPDataSinkChain() ) ;
	}

	void open_filtered_in_snp_data_sink() {
		reset_filtered_in_snp_data_sink() ;
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

	void open_filtered_out_snp_data_sink() {
		reset_filtered_out_snp_data_sink() ;
		if( m_options.check_if_option_was_supplied( "-write-snp-excl-list" ) && m_options.snp_excl_list_filename_mapper().output_filenames().size() > 0 ) {
			for( std::size_t i = 0; i < m_options.snp_excl_list_filename_mapper().output_filenames().size(); ++i ) {
				std::string const& filename = m_options.snp_excl_list_filename_mapper().output_filenames()[i] ;
				m_fltrd_out_snp_data_sink->add_sink( std::auto_ptr< genfile::SNPDataSink >( new SNPIDSink( filename ))) ;
			}
		}
		else {
			m_fltrd_out_snp_data_sink->add_sink( std::auto_ptr< genfile::SNPDataSink >( new genfile::TrivialSNPDataSink() )) ;
		}
	}

	void reset_filtered_out_snp_data_sink() {
		m_fltrd_out_snp_data_sink.reset( new genfile::SNPDataSinkChain() ) ;
	}

	void open_sample_row_sink() {
		m_fltrd_in_sample_sink.reset( new NullObjectSink< SampleRow >() ) ;
		if( m_options.output_sample_filename() != "" ) {
			m_backup_creator.backup_file_if_necessary( m_options.output_sample_filename() ) ;
			m_fltrd_in_sample_sink.reset( new SampleOutputFile< SimpleFileObjectSink< SampleRow > >( open_file_for_output( m_options.output_sample_filename() ))) ;
		}
		
		m_fltrd_out_sample_sink.reset( new NullObjectSink< SampleRow >() ) ;
		if( m_options.output_sample_excl_list_filename() != "" ) {
			m_fltrd_out_sample_sink.reset( new SampleIDSink( open_file_for_output( m_options.output_sample_excl_list_filename() ))) ;
		}
	}

	void open_snp_stats_sink( std::size_t index, GenRowStatistics const& snp_statistics ) {
		if( m_options.snp_stats_filename_mapper().output_filenames().size() == 0 ) {
			m_snp_stats_sink.reset( new statfile::TrivialBuiltInTypeStatSink() ) ;
		}
		else {
			assert( index < m_options.snp_stats_filename_mapper().output_filenames().size()) ;
			m_current_snp_stats_filename_index = index ;
			m_snp_stats_sink.reset( new statfile::RFormatStatSink( m_options.snp_stats_filename_mapper().output_filenames()[ index ] )) ;
		}
		m_snp_stats_sink->add_column( "" ) ;
		for( std::size_t i = 0; i < snp_statistics.size(); ++i ) {
			m_snp_stats_sink->add_column( snp_statistics.get_statistic_name( i )) ;
		}
	}

	void reset_snp_stats_sink() {
		m_snp_stats_sink.reset() ;
	}

	void open_sample_stats_sink() {
		if( m_options.output_sample_stats_filename() == "" ) {
			m_sample_stats_sink.reset( new statfile::TrivialBuiltInTypeStatSink() ) ;
		}
		else {
			m_sample_stats_sink.reset( new statfile::RFormatStatSink( m_options.output_sample_stats_filename() )) ;
		}
		m_sample_stats_sink->add_column( "" ) ;
		for( std::size_t i = 0; i < m_sample_statistics.size(); ++i ) {
			m_sample_stats_sink->add_column( m_sample_statistics.get_statistic_name( i )) ;
		}
	}
	
	void reset_sample_stats_sink() {
		m_sample_stats_sink.reset() ;
	}

	void construct_snp_statistics() {
		GenRowStatisticFactory::add_statistics( m_options.row_statistics_specs(), m_snp_statistics ) ;
	}

	void construct_sample_statistics() {
		std::vector< std::string > sample_statistics_specs = string_utils::split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-sample-stats-columns" ), "," ) ;
		SampleRowStatisticFactory::add_statistics( sample_statistics_specs, m_sample_statistics ) ;
	}

	void construct_snp_filter() {
		Timer timer ;
		m_logger << "(Constructing SNP filter...)" << std::flush ;
		
		std::auto_ptr< AndRowCondition > snp_filter( new AndRowCondition() ) ;

		if( m_options.check_if_option_was_supplied( "-hwe" ) ) {
			add_one_arg_condition_to_filter< StatisticLessThan >( *snp_filter, "HWE", m_options.get_value< double >( "-hwe" )) ;
		}

		if( m_options.check_if_option_was_supplied( "-info" ) ) {
			add_two_arg_condition_to_filter< StatisticInInclusiveRange >( *snp_filter, "information", m_options.get_values< double >( "-info" )) ;
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

		if( m_options.check_if_option_was_supplied( "-snpid-filter" ) ) {
			std::string expression = m_options.get_value< std::string >( "-snpid-filter" ) ;
			std::auto_ptr< RowCondition > snpid_condition( new SNPIDMatchesCondition( expression )) ;
			std::auto_ptr< RowCondition > real_snpid_condition( new NotRowCondition( snpid_condition )) ;
			snp_filter->add_subcondition( real_snpid_condition ) ;
		}

		if( m_options.check_if_option_was_supplied( "-snp-excl-list" ) ) {
			std::vector< std::string > filenames = m_options.get_values< std::string >( "-snp-excl-list" ) ;
			std::auto_ptr< RowCondition > snp_incl_condition( new SNPInListCondition( filenames )) ;
			std::auto_ptr< RowCondition > snp_excl_condition( new NotRowCondition( snp_incl_condition )) ;
			snp_filter->add_subcondition( snp_excl_condition ) ;
		}
		
		m_snp_filter = snp_filter ;
		m_snp_filter_failure_counts.resize( m_snp_filter->number_of_subconditions(), 0 ) ;
		
		m_logger << " (" << std::fixed << std::setprecision(1) << timer.elapsed() << "s)\n" ;
	}

	void construct_sample_filter() {
		Timer timer ;
		m_logger << "(Constructing sample filter...)" << std::flush  ;
		
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
		
		m_logger << " (" << std::fixed << std::setprecision(1) << timer.elapsed() << "s)\n" ;
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
	
	void load_sample_rows( std::size_t const expected_number_of_samples ) {
		
		try {
			unsafe_load_sample_rows( expected_number_of_samples ) ;
		}
		catch( ConditionValueNotFoundException const& ) {
			m_logger << "\n\n!! ERROR: The input sample file must contain entries for all values used to filter on.\n"
				<< "!! This includes \"missing\" and \"heterozygosity\".\n" ;
			throw ;
		}
		catch( genfile::MalformedInputError const& e ) {
			m_logger << "\n\n!! ERROR (" << e.what() << "): the sample file \"" << e.source() << " is malformed on line "
				<< e.line() ;
			if( e.has_column() ) {
				m_logger << ", column " << e.column() ;
			}
			m_logger << ".  Quitting.\n" ;
			throw ;
		}
	}
	
	void unsafe_load_sample_rows( std::size_t const expected_number_of_samples ) {
		if( m_options.input_sample_filename() != "" ) {
			genfile::CategoricalCohortIndividualSource sample_source( m_options.input_sample_filename() ) ;
			SampleRow sample_row ;
			for( std::size_t i = 0; i < sample_source.get_number_of_individuals(); ++i ) {
				sample_row.read_ith_sample_from_source( i, sample_source ) ;
				m_sample_rows.push_back( sample_row ) ;
				if( !m_sample_filter->check_if_satisfied( sample_row )) {
					m_indices_of_filtered_out_samples.push_back( m_sample_rows.size() - 1 ) ;
				}
			}
			if( m_sample_rows.size() != expected_number_of_samples ) {
				throw NumberOfSamplesMismatchException() ;
			}
		}
		else {
			m_sample_rows.resize( expected_number_of_samples ) ;
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
		if( m_options.check_if_option_was_supplied_in_group( "Sample filtering options") && m_options.output_sample_excl_list_filename() == "" && m_options.gen_filename_mapper().output_filenames().size() == 0 && m_options.snp_stats_filename_mapper().output_filenames().size() == 0 ) {
			m_errors.push_back( "You have specified sample filters, but not an output sample exclusion list, nor any output GEN files, nor any output SNP statistic files." ) ;
		}
		if( m_snp_filter->number_of_subconditions() > 0 && m_options.snp_excl_list_filename_mapper().output_filenames().size() == 0 && m_options.gen_filename_mapper().output_filenames().size() == 0 ) {
			m_errors.push_back( "You have specified SNP filters, but no output SNP exclusion list or output GEN files.\n" ) ;
		}
	}
	
	void check_for_warnings() {
		if( (m_options.output_sample_filename() != "" || m_options.output_sample_stats_filename() != "" || m_options.output_sample_excl_list_filename() != "" ) && m_options.gen_filename_mapper().input_files().size() != 22 ) {
			m_warnings.push_back( "You are outputting a sample, sample statistic, or sample exclusion file, but the number of gen files is not 22.\n"
			"   (I suspect there is not the whole genomes' worth of data?)" ) ;
		}
		if( m_options.output_sample_stats_filename() != "" && m_options.input_sample_filename() == "" ) {
			m_warnings.push_back( "You are outputting a sample statistic file, but no input sample file has been supplied.\n"
			"   Statistics will be output but the ID fields will be left blank.") ;
		}
		if( m_options.gen_filename_mapper().output_filenames().size() == 0 && m_options.output_sample_filename() == "" && m_options.snp_stats_filename_mapper().output_filenames().size() == 0 && m_options.output_sample_stats_filename() == "" && m_options.output_sample_excl_list_filename() == "" && m_options.snp_excl_list_filename_mapper().output_filenames().size() == 0 ) {
			m_warnings.push_back( "You have not specified any output files.  This will produce only logging output." ) ;
		}
		if( ((m_options.gen_filename_mapper().output_filenames().size() > 0) || ( m_options.snp_excl_list_filename_mapper().output_filenames().size() > 0)) &&  (m_snp_filter->number_of_subconditions() == 0) && (m_sample_filter->number_of_subconditions() == 0)) {
			m_warnings.push_back( "You have specified output GEN (or snp exclusion) files, but no filters.\n"
				" To convert GEN file formats, use gen-convert." ) ;
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
	std::auto_ptr< statfile::BuiltInTypeStatSink > m_snp_stats_sink ;
	std::auto_ptr< statfile::BuiltInTypeStatSink > m_sample_stats_sink ;

	std::vector< std::size_t > m_snp_filter_failure_counts ;
	std::vector< std::size_t > m_sample_filter_failure_counts ;
	
	bool m_ignore_warnings ; 
	
	std::size_t m_current_snp_stats_filename_index ;

	std::auto_ptr< AndRowCondition > m_snp_filter ;
	std::auto_ptr< AndRowCondition > m_sample_filter ;

	GenRowStatistics m_snp_statistics ;
	SampleRowStatistics m_sample_statistics ;
	
	std::vector< std::size_t > m_indices_of_filtered_out_samples ;
	
	ToNumberedFileBackupCreator m_backup_creator ;
	
	std::vector< std::string > m_warnings ;
	std::vector< std::string > m_errors ;
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
