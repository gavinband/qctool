/*
 * This program, gen-case-control-test
 */

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <numeric>
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
#include "parse_utils.hpp"
#include "Timer.hpp"

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
	std::vector< std::string > result( filenames ) ;
	if( filenames.size() != 2 ) {
		throw OptionValueInvalidException( option_name, filenames, "Two files, a sample and a gen file, must be specified for this option." ) ;
	}
	// swap the two options if it looks like the gen file is last.	
	if( !genfile::filename_indicates_gen_or_bgen_format( filenames[1] ) && genfile::filename_indicates_gen_or_bgen_format( filenames[0] )) {
		// Guess that sample file is listed second.
		std::swap( result[0], result[1] ) ;
	}
	return result ;
}

/*
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
		SampleRow const& m_sample_row ;
		double m_probability_of_case_status ;
} ;
*/
/*
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
		m_samples_by_probability.insert( std::make_pair( probability_sample_is_case, sample_index )) ;
	}

	std::vector< SampleSpec > m_samples ;
	std::multimap< double, std::size_t > m_samples_by_probability ;
} ;
*/

struct GenCaseControlProcessor
{
public:
	static void declare_options( OptionProcessor & options ) {
		
		// File options		
	    options[ "-cases" ]
	        .set_description( "Path of gen file to input" )
			.set_is_required()
			.set_takes_values()
			.set_maximum_number_of_repeats( 1 )
			.set_number_of_values_per_use( 2 )
			.add_value_preprocessor( &try_to_put_sample_file_first ) ;

	    options[ "-controls" ]
	        .set_description( "Path of gen file to input" )
			.set_is_required()
			.set_takes_values()
			.set_maximum_number_of_repeats( 1 )
			.set_number_of_values_per_use( 2 )
			.add_value_preprocessor( &try_to_put_sample_file_first ) ;

		options [ "--force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
	}

	GenCaseControlProcessor( OptionProcessor const& options )
		: m_cout( std::cout.rdbuf() ),
		  m_options( options )
	{
		write_start_banner( m_cout ) ;
		setup() ;
	}

	~GenCaseControlProcessor() {
		write_end_banner( m_cout ) ;
	}
	
private:
	void setup() {
		try {
			get_required_filenames() ;
			open_gen_files() ;
			open_sample_files() ;
			// m_case_control_test.reset( CaseControlTest::create( m_control_sample_source, m_control_gen_input_chain, m_case_sample_source, m_case_gen_input_chain )) ;
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

		control_filenames = m_options.get_values< std::string >( "-controls" ) ;
		assert( control_filenames.size() == 2 ) ;
		add_filename( &m_control_sample_filenames, control_filenames[0] ) ;
		expand_and_add_filename( &m_control_gen_filenames, control_filenames[1] ) ;

		case_filenames = m_options.get_values< std::string >( "-cases" ) ;
		assert( case_filenames.size() == 2 ) ;
		add_filename( &m_case_sample_filenames, case_filenames[0] ) ;
		expand_and_add_filename( &m_case_gen_filenames, case_filenames[1] ) ;
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
		open_gen_files( m_control_gen_filenames, m_control_gen_input_chain ) ;
		open_gen_files( m_case_gen_filenames, m_case_gen_input_chain ) ;

		if( timer.elapsed() > 1.0 ) {
			m_cout << "Opened " << m_control_gen_filenames.size() + m_case_gen_filenames.size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
		}
	}
	
	void open_gen_files( std::vector< std::string > const& filenames, std::auto_ptr< genfile::SNPDataSourceChain >& chain ) {
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
		// m_control_sample_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_control_sample_filenames[0] ))) ;
		// m_case_sample_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_case_sample_filenames[0] ))) ;
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
			oStream << "Control GEN files:\n" ;
			print_gen_files( oStream, m_control_gen_filenames, m_control_gen_input_chain ) ;
			oStream << "Case GEN files:\n" ;
			print_gen_files( oStream, m_case_gen_filenames, m_case_gen_input_chain ) ;
			oStream << std::string( 72, '=' ) << "\n\n" ;
		}
		catch (...) {
			oStream << std::string( 72, '=' ) << "\n\n" ;
			throw ;
		}
	}

	void print_gen_files( std::ostream& oStream, std::vector<std::string> const& filenames, std::auto_ptr< genfile::SNPDataSourceChain > const& chain ) const {
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			oStream
				<< "  (" << std::setw(6) << chain->number_of_snps_in_source( i ) << " snps)  "
				<< "\"" << std::setw(20) << filenames[i] << "\"\n" ;
		}
		oStream << "  (total " << chain->total_number_of_snps() << " snps).\n\n" ;
	}

	void do_checks() {
	}
	
	void process() {
		unsafe_process() ;
	}

private:

	struct SNPsHaveSamePositionChecker
	{
		SNPsHaveSamePositionChecker( GenRow const& row ): m_row( row ) {}
		bool operator()( std::string const&, std::string const&, int SNP_position, char, char ) const {
			return m_row.SNP_position() == SNP_position ;
		}
	private:
		GenRow const& m_row ;
	} ;

	void unsafe_process() {
		Timer timer ;

		GenRow control_row, case_row ;
		control_row.set_number_of_samples( m_control_gen_input_chain->number_of_samples() ) ;
		case_row.set_number_of_samples( m_case_gen_input_chain->number_of_samples() ) ;

		m_cout << "Processing SNPs...\n" ;

		double last_time = -5.0 ;
		std::size_t number_of_case_snps_processed = 0 ;
		std::size_t number_of_control_snps_matched = 0 ;
		while( read_snp( *m_case_gen_input_chain, case_row )) {
			++number_of_case_snps_processed ;
			std::size_t number_of_matching_snps = 0 ;
			if( read_next_matching_snp( *m_control_gen_input_chain, control_row, SNPsHaveSamePositionChecker( case_row ))) {
				++number_of_matching_snps ;
				++number_of_control_snps_matched ;
				// I think we'd expect all the SNP identifying data to be the same for the two snps.
				// If there are mismatches, print out some information about them.
				compare_snps( case_row, control_row ) ;
				
				// print a progress message every second.
				double time_now = timer.elapsed() ;
				if( 1 ) {//(time_now - last_time >= 1.0) || (number_of_case_snps_processed == m_case_gen_input_chain->total_number_of_snps()) ) {
					print_progress( time_now ) ;
					last_time = time_now ;
				}
			}
			if( number_of_matching_snps != 1 ) {
				m_cout << "\nA strange number (" << number_of_matching_snps << ") of control snps matched the case snp "
					<< case_row.SNPID() << " "
					<< case_row.RSID() << " "
					<< case_row.SNP_position()
					<< ".\n" ;
			}
		}

		std::cerr << "\nProcessed case / control data (" << number_of_case_snps_processed << " case SNPs, " << number_of_control_snps_matched << " matched control SNPs) in " << timer.elapsed() << " seconds.\n" ;
	
		m_cout << "Post-processing..." << std::flush ;
		timer.restart() ;
		// Close the output gen file(s) now
		close_all_files() ;
		std::cerr << " (" << timer.elapsed() << "s)\n" ;
	}
	
	void compare_snps( GenRow const& case_row, GenRow const& control_row ) {
		if( case_row.SNPID() != control_row.SNPID() ) {
			m_cout << "\nMatching rows have differing SNPIDs " << case_row.SNPID() << " and " << control_row.SNPID() << ".\n" ;
		}
		if( case_row.RSID() != control_row.RSID() ) {
			m_cout << "\nMatching rows have differing RSIDs " << case_row.RSID() << " and " << control_row.RSID() << ".\n" ;
		}
		if( case_row.SNP_position() != control_row.SNP_position() ) {
			m_cout << "\nMatching rows have different positions " << case_row.SNP_position() << " and " << control_row.SNP_position() << ".\n" ;
		}
		if( case_row.first_allele() != control_row.first_allele() || case_row.second_allele() != control_row.second_allele() ) {
			m_cout 	<< "\nMatching rows have differing alleles "
					<< case_row.first_allele() << " " << case_row.second_allele()
					<< " and " << control_row.first_allele() << " " << control_row.second_allele() << ".\n" ;
		}
	}

	void print_progress( double time_now ) {
		double case_progress = (static_cast< double >( m_case_gen_input_chain->number_of_snps_read() ) / m_case_gen_input_chain->total_number_of_snps()) ;
		double control_progress = (static_cast< double >( m_control_gen_input_chain->number_of_snps_read() ) / m_control_gen_input_chain->total_number_of_snps()) ;
		m_cout
			<< "\r"
			<< get_progress_bar( 30, case_progress )
			<< " (" << m_case_gen_input_chain->number_of_snps_read() << "/" << m_case_gen_input_chain->total_number_of_snps() << ")"
			<< " "
			<< get_progress_bar( 30, control_progress )
			<< " (" << m_control_gen_input_chain->number_of_snps_read() << "/" << m_control_gen_input_chain->total_number_of_snps() << ")"
			<< " (" << std::fixed << std::setprecision(1) << time_now << "s)"
			<< std::string( std::size_t(5), ' ' )
			<< std::flush ;
	}

	std::string get_progress_bar( std::size_t width, double progress ) {
		progress = std::min( std::max( progress, 0.0 ), 1.0 ) ;
		std::size_t visible_progress = progress * width ;
		return
			"["
			+ std::string( std::size_t( visible_progress ), '*' )
			+ std::string( std::size_t( width - visible_progress ), ' ' )
			+ "]" ;
	}

	typedef boost::function< bool( std::string const&, std::string const&, uint32_t, char, char ) > SNPMatcher ;

	bool read_snp( genfile::SNPDataSource& snp_data_source, GenRow& row ) {
		return snp_data_source.read_snp(
			boost::bind< void >( &GenRow::set_number_of_samples, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNPID, &row, _1 ),
			boost::bind< void >( &GenRow::set_RSID, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNP_position, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele1, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele2, &row, _1 ),
			boost::bind< void >( &GenRow::set_genotype_probabilities, &row, _1, _2, _3, _4 )
		) ;
	}

	bool read_next_matching_snp( genfile::SNPDataSource& snp_data_source, GenRow& row, SNPMatcher const& snp_matcher ) {
		return snp_data_source.read_next_matching_snp(
			boost::bind< void >( &GenRow::set_number_of_samples, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNPID, &row, _1 ),
			boost::bind< void >( &GenRow::set_RSID, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNP_position, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele1, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele2, &row, _1 ),
			boost::bind< void >( &GenRow::set_genotype_probabilities, &row, _1, _2, _3, _4 ),
			snp_matcher
		) ;
	}

	void close_all_files() {
		m_case_gen_input_chain.reset() ;
		m_control_gen_input_chain.reset() ;
	}
	
private:
	
	std::ostream m_cout ;
	
	std::auto_ptr< genfile::SNPDataSourceChain > m_case_gen_input_chain ;
	std::auto_ptr< genfile::SNPDataSourceChain > m_control_gen_input_chain ;

	OptionProcessor const& m_options ;
	
	std::vector< std::string > m_case_gen_filenames ;
	std::vector< std::string > m_case_sample_filenames ;
	std::vector< std::string > m_control_gen_filenames ;
	std::vector< std::string > m_control_sample_filenames ;

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
	}
	catch( std::exception const& e )
	{
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}
