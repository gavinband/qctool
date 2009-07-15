/*
 * This program, gen-compare, compares two gen files (or numbered collections of gen files)
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
#include "string_utils.hpp"

std::vector< std::string > expand_filename_wildcards( std::string const& option_name, std::vector< std::string > const& filenames ) ;
void check_files_are_readable( std::string const& option_name, std::vector< std::string > const& filenames ) ;

struct GenCompareProcessorException: public GToolException
{
	char const* what() const throw() { return "GenCompareProcessorException" ; }
} ;

struct GenCompareUsageError: public GenCompareProcessorException
{
	char const* what() const throw() { return "GenCompareUsageError" ; }
} ;

// thrown to indicate that the numbers of input and output files specified on the command line differ.
struct GenCompareFileCountMismatchError: public GenCompareProcessorException
{
	char const* what() const throw() { return "GenCompareFileCountMismatchError" ; }
} ;

// thrown to indicate that an output file with wildcard appeared, but the corresponding input
// file had no wildcard.
struct GenCompareFileWildcardMismatchError: public GenCompareProcessorException
{
	char const* what() const throw() { return "GenCompareFileWildcardMismatchError" ; }
} ;

struct GenCompareProcessor
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

	    options[ "--og" ]
	        .set_description( "Path of gen file to output" )
			.set_takes_values()
			.set_maximum_number_of_repeats( 100 ) ;
        	.add_value_checker( &check_files_dont_exist ) ;
	}

	GenCompareProcessor( OptionProcessor const& options )
		: m_cout( std::cout.rdbuf() ),
		  m_options( options )
	{
		setup() ;
	}
	
private:
	void setup() {
		get_required_filenames() ;
		try {
			open_gen_input_files() ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_cout << "The GEN files specified did not all have the same sample size.\n" ;
			throw ;
		}
	}

	void get_required_filenames() {
		std::vector< std::string > gen_filenames, gen_output_filenames ;
		if( m_options.check_if_option_was_supplied( "--g" ) ) {
			gen_filenames = m_options.get_values< std::string >( "--g" ) ;
		}
		if( m_options.check_if_option_was_supplied( "--og" ) ) {
			gen_output_filenames = m_options.get_value< std::string >( "--og" ) ;
		}

		if( gen_filenames.size() != gen_output_filenames.size() ) {
			throw GenCompareFileCountMismatchError() ;
		}

		add_expanded_filenames( gen_filenames, gen_output_filenames ) ;

		// For each entry in m_gen_filenames, there should now be a corresponding (index, output filename)
		// entry in m_gen_output_filenames.
		assert( m_gen_filenames.size() == m_gen_output_filenames.size() ) ;
	}

	void add_expanded_filenames( std::vector< std::string > const& gen_filenames, std::vector< std::string > const& gen_output_filenames ) {		
		for( std::size_t i = 0; i < gen_filenames.size(); ++i ) {
			add_expanded_filename( gen_filenames[i], gen_output_filenames[i] ) ;
		}
	}

	void add_expanded_filename( std::string const& input_filename, std::string const& output_filename ) {
		bool input_file_has_wildcard = ( input_filename.find( '#' ) != std::string::npos ) ;
		bool output_file_has_wildcard = ( output_filename.find( '#' ) != std::string::npos ) ;
		
		if( output_file_has_wildcard && !input_file_has_wildcard ) {
			throw GenCompareFileWildcardMismatchError() ;
		}

		if( input_file_has_wildcard ) {
			std::pair< std::vector< std::string >, std::vector< std::string > >
				expanded_input_filename = find_files_matching_path_with_wildcard( input_filename, '#' ) ;

			// we only use matching filenames if the match is a number from 1 to 100
			// For such filenames, we place a corresponding filename in the list of output files.
			for( std::size_t j = 0; j < expanded_input_filename.first.size(); ) {
				if( check_if_string_is_a_number_from_1_to_100( expanded_input_filename.second[j] )) {
					add_input_and_corresponding_output_filename( expanded_input_filename.first[j], output_filename, expanded_input_filename.second[j] ) ;
				} else {
					++j ;
				}
			}
		}
		else {
			add_input_and_corresponding_output_filename( expanded_input_filename.first[j], output_filename, "" ) ;
		}
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

	void add_input_and_corresponding_output_filename( std::string const& input_filename, std::string const& output_filename, std::string const& input_filename_wildcard_match ) {
		std::string expanded_output_filename = output_filename ;
		std::string::iterator wildcard_pos = output_filename.find( '#' ) ;
		if( wildcard_pos != std::string::npos ) {
			expanded_output_filename.replace( wildcard_pos, wildcard_pos + 1, input_filename_wildcard_match ) ;			
		}
		m_gen_input_filenames.push_back( input_filename ) ;
		m_gen_output_filenames[ m_gen_input_filenames.size() ] = expanded_output_filename ;
	}

	void open_gen_input_files() {
		#ifdef HAVE_BOOST_TIMER
			boost::timer timer ;
		#endif

		m_input_chain.reset( new genfile::SNPDataSourceChain() ) ;
		m_input_chain.set_moved_to_next_source_callback( boost::bind( &GenCompareProcessor::move_to_next_output_file, this )) ;
		
		std::size_t number_of_snps = 0 ;
		for( std::size_t i = 0; i < m_gen_filenames.size(); ++i ) {
			add_gen_file_to_chain( *chain, m_gen_filenames[i] ) ;
			m_gen_file_snp_counts.push_back( chain->total_number_of_snps() - number_of_snps ) ;
			number_of_snps = chain->total_number_of_snps ;
		}
		
		#ifdef HAVE_BOOST_TIMER
			if( timer.elapsed() > 1.0 ) {
				m_cout << "Opened " << m_gen_filenames.size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
			}
		#endif
	}

	void add_gen_file_to_chain( genfile::SNPDataSourceChain& chain, std::string const& filename ) const {
		#ifdef HAVE_BOOST_TIMER
			boost::timer file_timer ;
		#endif
			m_cout << "(Opening gen file \"" << filename << "\"...)" << std::flush ;
			try {
				std::auto_ptr< genfile::SNPDataSource > snp_data_source( genfile::SNPDataSource::create( filename )) ;
				chain->add_source( snp_data_source ) ;
			}
			catch ( genfile::FileHasTwoConsecutiveNewlinesError const& e ) {
				std::cerr << "\n!!ERROR: a GEN file was specified having two consecutive newlines.\n"
					<< "!! NOTE: popular editors, such as vim and nano, automatically add an extra newline to the file (which you can't see).\n"
					<< "!!     : Please check that each SNP in the file is terminated by a single newline.\n" ;
				throw ;
			}
		#ifdef HAVE_BOOST_TIMER
			m_cout << " (" << file_timer.elapsed() << "s\n" ;
		#else
			m_cout << "\n" ;
		#endif			
	}

	void open_gen_output_files() {
		m_output_chain.reset( new genfile::SNPDataSinkChain() ) ;
		output_filenames_t::const_iterator
			i = m_gen_output_filenames.begin(),
			end_i = m_gen_output_filenames.end() ;
		for( ; i != end_i; ++i ) {
			add_gen_file_to_chain( *chain, i->second ) ;
		}
	}

	void add_gen_file_to_chain( genfile::SNPDataSinkChain& chain, std::string const& filename ) const {
			m_cout << "(Opening output gen file \"" << filename << "\"...)" << std::flush ;
			std::auto_ptr< genfile::SNPDataSource > snp_data_sink( genfile::SNPDataSink::create( filename )) ;
			chain->add_sink( snp_data_sink ) ;
	}

	void move_to_next_output_file() {
		m_output_chain.move_to_next_sink() ;
	}


public:
	
	void write_banner( std::ostream& oStream ) const {
		oStream << "\nWelcome to gen-compare\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}
	
	void write_preamble( std::ostream& oStream ) const {
		oStream << std::string( 72, '=' ) << "\n\n" ;
		try {
			oStream << std::setw(30) << "Input and output GEN files:" ;
			for( std::size_t i = 0; i < m_gen_filenames.size(); ++i ) {
				if( i > 0 ) {
					oStream << std::string( 30, ' ' ) ;
				}
				oStream
					<< "  (" << std::setw(6) << m_gen_file_snp_counts[i] << " snps)  "
					<< "\"" << std::setw(20) << m_gen_filenames[i] << "\""
					<< " -> \"" << m_gen_output_filenames[i] << "\"\n";
			}
			oStream << std::string( 30, ' ' ) << "  (total " << m_total_number_of_snps << " snps).\n" ;

			if( !m_errors.empty() ) {
				for( std::size_t i = 0; i < m_errors.size(); ++i ) {
					oStream << "!! ERROR: " << m_errors[i] << "\n\n" ;
				}
				oStream << "!! Please correct the above errors and re-run gen-compare.\n\n" ;
				throw GenCompareUsageError ;
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
			output_filenames_t::const_iterator
				j = m_gen_output_filenames.begin(),
				end_j = m_gen_output_filenames.end() ;
			for( ; j != end_j; ++j ) {
				if( strings_are_nonempty_and_equal( m_gen_filenames[i], j->second )) {
					m_errors.push_back( "Output GEN file \"" + j->second +"\" also specified as input GEN file." ) ;
					break ;
				}
			}
		}
	}
	
	bool strings_are_nonempty_and_equal( std::string const& left, std::string const& right ) {
		return (!left.empty()) && (!right.empty()) && (left == right) ;
	}

	
	void process() {
		unsafe_process() ;
	}

private:

	void unsafe_process() {
		process_gen_rows() ;
	}

	void process_gen_rows() {
	#ifdef HAVE_BOOST_TIMER
		boost::timer timer ;
	#endif

		open_gen_output_files() ;

		GenRow row ;
		row.set_number_of_samples( m_number_of_samples ) ;

		while( m_input_chain.read_snp(
			boost::bind< void >( &GenRow::set_number_of_samples, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNPID, &row, _1 ),
			boost::bind< void >( &GenRow::set_RSID, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNP_position, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele1, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele2, &row, _1 ),
			boost::bind< void >( &GenRow::set_genotype_probabilities, &row, _1, _2, _3, _4 )
		)) {
			process_gen_row( row, ++m_total_number_of_snps ) ;
		}

	#ifdef HAVE_BOOST_TIMER
		std::cerr << "gen-select: processed GEN file(s) (" << m_total_number_of_snps << " rows) in " << timer.elapsed() << " seconds.\n" ;
	#endif
	
		m_cout << "Post-processing (updating file header, compression)..." << std::flush ;
		timer.restart() ;
		// Close the output gen file(s) now
		close_gen_output_files() ;
	#ifdef HAVE_BOOST_TIMER
		std::cerr << " (" << timer.elapsed() << "s)\n" ;
	#else
		std::cerr << "\n" ;
	#endif
	}

	void process_gen_row( GenRow const& row, std::size_t row_number ) {
		output_gen_row( row ) ;
	}

	void output_gen_row( GenRow const& row ) {
		m_output_chain->write_snp(
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

private:
	
	std::ostream m_cout ;
	
	std::auto_ptr< genfile::SNPDataSourceChain > m_input_chain ;
	std::auto_ptr< genfile::SNPDataSinkChain > m_output_chain ;
	
	OptionProcessor const& m_options ;
	
	std::vector< std::size_t > m_gen_file_snp_counts ;
	std::size_t m_total_number_of_snps ;
	std::size_t m_number_of_samples_from_gen_file ;
	
	std::vector< std::string > m_gen_filenames ;
	typedef std::map< std::size_t, std::string > output_filenames_t ;
	output_filenames_t m_gen_output_filenames ;

	std::vector< std::string > m_errors ;
} ;


int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		GenCompareProcessor::declare_options( options ) ;
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
		GenCompareProcessor processor( options ) ;
		processor.write_banner( std::cout ) ;
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



