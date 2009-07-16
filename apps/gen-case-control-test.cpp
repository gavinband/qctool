/*
 * This program, gen-convert, compares two gen files (or numbered collections of gen files)
 * on a SNP-by-SNP basis.
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
#include "Whitespace.hpp"
#include "FileUtil.hpp"
#include "GenRowSource.hpp"
#include "GenRowSink.hpp"
#include "SNPDataSource.hpp"
#include "SNPDataSourceChain.hpp"
#include "SNPDataSink.hpp"
#include "SNPDataSinkChain.hpp"
#include "string_utils.hpp"
#include "wildcard.hpp"

namespace globals {
	std::string const program_name = "gen-case-control-test" ;
}

struct GenCaseControlProcessorException: public GToolException
{
	char const* what() const throw() { return "GenCaseControlProcessorException" ; }
} ;

struct GenCaseControlUsageError: public GenCaseControlProcessorException
{
	char const* what() const throw() { return "GenCaseControlUsageError" ; }
} ;

// thrown to indicate that the numbers of files specified for -cases or -controls
// is not two
struct GenCaseControlFileCountError: public GenCaseControlProcessorException
{
	char const* what() const throw() { return "GenCaseControlFileCountError" ; }
} ;

// thrown to indicate that an output file with wildcard appeared, but the corresponding input
// file had no wildcard.
struct GenCaseControlFileWildcardMismatchError: public GenCaseControlProcessorException
{
	char const* what() const throw() { return "GenCaseControlFileWildcardMismatchError" ; }
} ;

std::vector< std::string > try_to_put_sample_file_first( std::string const& option_name, std::vector< std::string > const& filenames ) {
	if( filenames.size() != 2 ) {
		throw OptionValueInvalidException( "Two files, a sample and a gen file, must be specified for this option." ) ;
	}
	// swap the two options if it looks like the gen file is last.	
	if( !genfile::filename_indicates_gen_or_bgen_format( filenames[1] ) && genfile::filename_indicates_gen_or_bgen_format( filenames[0] )) {
		// Guess that sample file is listed second.
		std::swap( filenames.begin(), filenames.begin() + 1 ) ;
	}
}


struct SampleSpec
{
	public:
		SampleSpec( SampleRow const& row, double probability_of_case_status )
			: m_sample_row( row ), m_probability_of_case_status( probability_of_case_status )
		{
		}
		
		SampleRow const& sample_row() const { return m_sample_row ; }
		double probability_of_case_status() const { return m_probability_of_case_status ; }

	private:
		SampleRow m_sample_row ;
		double m_probability_of_case_status ;
} ;


struct CaseControlTest
{
public:
	static std::auto_ptr< CaseControlTest > create(
		std::auto_ptr< ObjectSource< SampleRow > > case_sample_source,
		std::auto_ptr< genfile::SNPDataSource > control_snp_data_source,
		std::auto_ptr< ObjectSource< SampleRow > > control_sample_source,
		std::auto_ptr< genfile::SNPDataSource > case_snp_data_source
	) {
	}

private:
	
	void add_sample( SampleRow const& row, double probability_sample_is_case ) {
		std::size_t sample_index = m_samples.size() ;
		m_samples.push_back( SampleSpec( row, probability_sample_is_case )) ;
		m_samples_by_probability.insert( std::make_pair( probability_sample_is_case, sample_index ))
	}

	std::vector< SampleSpec > m_samples ;
	std::multimap< double, std::size_t > m_samples_by_probability ;
} ;


struct GenCaseControlProcessor
{
public:
	static void declare_options( OptionProcessor & options ) {
		
		// File options		
	    options[ "-cases" ]
	        .set_description( "Path of gen file to input" )
			.set_required()
			.set_takes_values()
			.set_maximum_number_of_repeats( 1 )
			.set_number_of_values_per_use( 2 )
			.set_option_preprocessor( &try_to_put_sample_file_first ) ;

	    options[ "-controls" ]
	        .set_description( "Path of gen file to input" )
			.set_required()
			.set_takes_values()
			.set_maximum_number_of_repeats( 1 )
			.set_number_of_values_per_use( 2 )
			.set_option_preprocessor( &try_to_put_sample_file_first ) ;

		options [ "--force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
	}

	GenCaseControlProcessor( OptionProcessor const& options )
		: m_cout( std::cout.rdbuf() ),
		  m_options( options )
	{
		setup() ;
	}
	
private:
	void setup() {
		write_banner( m_cout ) ;
		try {
			get_required_filenames() ;
			open_gen_files() ;
			open_sample_files() ;
			m_case_control_test.reset( CaseControlTest::create( m_control_sample_source, m_control_gen_input_chain, m_case_sample_source, m_case_gen_input_chain )) ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_cout << "The GEN files specified did not all have the same sample size.\n" ;
			throw ;
		}
		catch ( genfile::FileHasTwoConsecutiveNewlinesError const& e ) {
			m_cout << "\n!!ERROR: a GEN file was specified having two consecutive newlines.\n"
				<< "!! NOTE: popular editors, such as vim and nano, automatically add an extra newline to the file (which you can't see).\n"
				<< "!!     : Please check that each SNP in the file is terminated by a single newline.\n" ;
			throw ;
		}
	}

	void get_required_filenames() {
		std::vector< std::string > case_filenames, control_filenames ;

		case_filenames = m_options.get_values< std::string >( "-cases" ) ;
		expand_and_add_filename( &m_control_gen_filenames, control_filenames[0] ) ;
		add_filename( &m_control_sample_filenames, control_filenames[1] ) ;

		control_filenames = m_options.get_values< std::string >( "-controls" ) ;
		expand_and_add_filename( &m_case_gen_filenames, case_filenames[0] ) ;
		add_filename( &m_case_sample_filenames, case_filenames[1] ) ;
	}

	void expand_and_add_filename( std::vector< std::string >* filename_list_ptr, std::string const& filename ) {
		bool input_file_has_wildcard = ( filename.find( '#' ) != std::string::npos ) ;
		if( input_file_has_wildcard ) {
			std::pair< std::vector< std::string >, std::vector< std::string > >
				expanded_filename = find_files_matching_path_with_wildcard( filename, '#' ) ;

			// we only use matching filenames if the match is a number from 1 to 100
			// For such filenames, we add the filename to our list for cases.
			for( std::size_t j = 0; j < expanded_filename.first.size(); ++j ) {
				if( check_if_string_is_a_number_from_1_to_100( expanded_filename.second[j] )) {
					add_filename( filename_list_ptr, expanded_filename.first[j] ) ;
				}
			}
		}
		else {
			add_filename( filename_list_ptr, filename ) ;
		}
	}

	bool check_if_string_is_a_number_from_1_to_100( std::string const& a_string ) {
		return parse_integer_in_half_open_range( a_string, 1, 101 ) != 101 ;
	}

	void add_filename( std::vector< std::string >* filename_list_ptr, std::string const& filename ) {
		filename_list_ptr->push_back( filename ) ;
	}

	void open_gen_files() {
		Timer timer ;
		open_gen_files( m_control_gen_filenames, m_control_gen_input_chain )
		open_gen_files( m_case_gen_filenames, m_case_gen_input_chain )

		if( timer.elapsed() > 1.0 ) {
			m_cout << "Opened " << m_control_gen_filenames.size() + m_case_gen_filenames.size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
		}
	}
	
	void open_gen_files( std::vector< std::string > const& filenames, std::auto_ptr< SNPDataSourceChain >& chain ) {
		chain.reset( new genfile::SNPDataSourceChain() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			add_gen_file_to_chain( *chain, filenames[i] ) ;
		}
	}

	void add_gen_file_to_chain( genfile::SNPDataSourceChain& chain, std::string const& filename ) {
		Timer timer ;
		m_cout << "(Opening gen file \"" << filename << "\"...)" << std::flush ;
		std::auto_ptr< genfile::SNPDataSource > snp_data_source( genfile::SNPDataSource::create( filename )) ;
		chain.add_source( snp_data_source ) ;
		m_cout << " (" << timer.elapsed() << "s)\n" ;
	}

	void open_sample_files() {
		assert( m_control_sample_filenames.size() == 1 ) ;
		Timer timer ;
		m_control_sample_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_control_sample_filenames[0] ))) ;
		m_case_sample_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_case_sample_filenames[0] ))) ;
	}

public:
	
	void write_start_banner( std::ostream& oStream ) const {
		oStream << "\nWelcome to " << globals::program_name << "\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}
	
	void write_end_banner( std::ostream& oStream ) const {
		oStream << "\n"
			<< "Thank you for using " << globals::program_name << ".\n" ;
	}

	void write_preamble( std::ostream& oStream ) const {
		oStream << std::string( 72, '=' ) << "\n\n" ;
		try {
			oStream << std::setw(30) << "Input and output GEN files:" ;
			for( std::size_t i = 0; i < m_gen_input_filenames.size(); ++i ) {
				if( i > 0 ) {
					oStream << std::string( 30, ' ' ) ;
				}
				oStream
					<< "  (" << std::setw(6) << m_gen_file_snp_counts[i] << " snps)  "
					<< "\"" << std::setw(20) << m_gen_input_filenames[i] << "\"\n" ;
			}
			oStream << std::string( 30, ' ' ) << "  (total " << m_input_chain->total_number_of_snps() << " snps).\n\n" ;
			oStream << std::string( 72, '=' ) << "\n\n" ;
		}
		catch (...) {
			oStream << std::string( 72, '=' ) << "\n\n" ;
			throw ;
		}
	}

	void do_checks() {
	}
	
	void process() {
		unsafe_process() ;
	}

private:

	void unsafe_process() {
		process_gen_rows() ;
	}

	void process_gen_rows() {
		Timer timer ;

		GenRow row ;
		row.set_number_of_samples( m_number_of_samples_from_gen_file ) ;

		m_cout << "Converting GEN files...\n" ;

		double last_time = -5.0 ;
		std::size_t number_of_snps_processed = 0 ;
		while( read_snp(row) ) {
			++number_of_snps_processed ;
			// print a message every 5 seconds.

			double time_now = timer.elapsed() ;
			if( time_now - last_time >= 1.0 || number_of_snps_processed == m_input_chain->total_number_of_snps() ) {
				std::size_t progress = (static_cast< double >( number_of_snps_processed ) / m_input_chain->total_number_of_snps()) * 30.0 ;
				m_cout
					<< "\r["
					<< std::string( progress, '*' )
					<< std::string( 30 - progress, ' ' )
					<< "]"
					<< " ("
					<< number_of_snps_processed
					<< "/"
					<< m_input_chain->total_number_of_snps()
					<< " SNPs, "
					<< std::fixed << std::setprecision(1) << timer.elapsed()
					<< "s, "
					<< std::setw( 5 ) << std::fixed << std::setprecision(1) << (number_of_snps_processed / timer.elapsed())
					<< " SNPs/s)" 
					<< std::string( std::size_t(5), ' ' )
 					<< std::flush ;
				last_time = time_now ;
			}

			write_snp(row) ;
		}

		std::cerr << "Converted GEN file(s) (" << number_of_snps_processed << " SNPs) in " << timer.elapsed() << " seconds.\n" ;
	
		m_cout << "Post-processing (updating file header, compression)..." << std::flush ;
		timer.restart() ;
		// Close the output gen file(s) now
		close_all_files() ;
		std::cerr << " (" << timer.elapsed() << "s)\n" ;
	#else
		std::cerr << "\n" ;
	#endif
	}

	bool read_snp( GenRow& row ) {
		return m_input_chain->read_snp(
			boost::bind< void >( &GenRow::set_number_of_samples, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNPID, &row, _1 ),
			boost::bind< void >( &GenRow::set_RSID, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNP_position, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele1, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele2, &row, _1 ),
			boost::bind< void >( &GenRow::set_genotype_probabilities, &row, _1, _2, _3, _4 )
		) ;
	}

	bool write_snp( GenRow& row ) {
		return m_output_chain->write_snp(
			row.number_of_samples(),
			row.SNPID(),
			row.RSID(),
			row.SNP_position(),
			row.first_allele(),
			row.second_allele(),
			boost::bind< double >( &GenRow::get_AA_probability, &row, _1 ),
			boost::bind< double >( &GenRow::get_AB_probability, &row, _1 ),
			boost::bind< double >( &GenRow::get_BB_probability, &row, _1 )
		) ;
	}

	void close_all_files() {
		m_input_chain.reset() ;
		m_output_chain.reset() ;
	}
	
private:
	
	std::ostream m_cout ;
	
	std::auto_ptr< genfile::SNPDataSourceChain > m_input_chain ;
	std::auto_ptr< genfile::SNPDataSinkChain > m_output_chain ;
	std::string m_current_output_filename ;

	OptionProcessor const& m_options ;
	
	std::vector< std::size_t > m_gen_file_snp_counts ;
	std::size_t m_number_of_samples_from_gen_file ;
	
	std::vector< std::string > m_gen_input_filenames ;
	typedef std::vector< std::string > output_filenames_t ;
	output_filenames_t m_gen_output_filenames ;

	std::size_t m_output_file_index ;

	std::vector< std::string > m_errors ;
} ;


int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		GenCaseControlProcessor::declare_options( options ) ;
		options.process( argc, argv ) ;
    }
    catch( std::exception const& exception ) {
        std::cerr << "!! Error: " << exception.what() << ".\n";
        std::cerr << "Usage: " << globals::program_name << " [options]\n"
                << options
                << "\n" ;
        return -1 ;
    }

	// main section

	try	{
		GenCaseControlProcessor processor( options ) ;
		processor.do_checks() ;
		processor.write_preamble( std::cout ) ;
		processor.process() ;
		processor.write_postamble( std::cout ) ;
	}
	catch( std::exception const& e )
	{
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}
