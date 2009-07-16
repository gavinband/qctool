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
#include <boost/bind.hpp>
#include "GenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "OptionProcessor.hpp"
#include "FileUtil.hpp"
#include "SNPDataSource.hpp"
#include "SNPDataSourceChain.hpp"
#include "SNPDataSink.hpp"
#include "SNPDataSinkChain.hpp"
#include "string_utils.hpp"
#include "wildcard.hpp"

namespace globals {
	std::string const program_name = "gen-convert" ;
}

struct GenConvertProcessorException: public GToolException
{
	char const* what() const throw() { return "GenConvertProcessorException" ; }
} ;

struct GenConvertUsageError: public GenConvertProcessorException
{
	char const* what() const throw() { return "GenConvertUsageError" ; }
} ;

// thrown to indicate that the numbers of input and output files specified on the command line differ.
struct GenConvertFileCountMismatchError: public GenConvertProcessorException
{
	char const* what() const throw() { return "GenConvertFileCountMismatchError" ; }
} ;

// thrown to indicate that an output file with wildcard appeared, but the corresponding input
// file had no wildcard.
struct GenConvertFileWildcardMismatchError: public GenConvertProcessorException
{
	char const* what() const throw() { return "GenConvertFileWildcardMismatchError" ; }
} ;

struct GenConvertProcessor
{
public:
	static void declare_options( OptionProcessor & options ) {
		
		// File options		
	    options[ "--g" ]
	        .set_description( "Path of gen file to input" )
	        .set_is_required()
			.set_takes_values()
			.set_maximum_number_of_repeats( 100 ) ;

	    options[ "--og" ]
	        .set_description( "Path of gen file to output" )
			.set_takes_values()
			.set_maximum_number_of_repeats( 100 ) ;
			
		options [ "--force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
	}

	GenConvertProcessor( OptionProcessor const& options )
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
			open_gen_input_files() ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_cout << "The GEN files specified did not all have the same sample size.\n" ;
			throw ;
		}
		catch( GenConvertFileCountMismatchError const& ) {
			m_cout << "The number of files specified for --g must match that for --og.\n" ;
			throw ;
		}
		catch( GenConvertFileWildcardMismatchError const& ) {
			m_cout << "The output file contained a wildcard character, so the input should also.\n" ;
			throw ;
		}
	}

	void get_required_filenames() {
		std::vector< std::string > gen_filenames, gen_output_filenames ;
		if( m_options.check_if_option_was_supplied( "--g" ) ) {
			gen_filenames = m_options.get_values< std::string >( "--g" ) ;
		}
		if( m_options.check_if_option_was_supplied( "--og" ) ) {
			gen_output_filenames = m_options.get_values< std::string >( "--og" ) ;
		}

		if( gen_filenames.size() != gen_output_filenames.size() ) {
			throw GenConvertFileCountMismatchError() ;
		}

		add_expanded_filenames( gen_filenames, gen_output_filenames ) ;

		// For each entry in m_gen_input_filenames, there should now be a corresponding (index, output filename)
		// entry in m_gen_output_filenames.
		assert( m_gen_input_filenames.size() == m_gen_output_filenames.size() ) ;
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
			throw GenConvertFileWildcardMismatchError() ;
		}

		if( input_file_has_wildcard ) {
			std::pair< std::vector< std::string >, std::vector< std::string > >
				expanded_input_filename = find_files_matching_path_with_wildcard( input_filename, '#' ) ;

			// we only use matching filenames if the match is a number from 1 to 100
			// For such filenames, we place a corresponding filename in the list of output files.
			for( std::size_t j = 0; j < expanded_input_filename.first.size(); ++j ) {
				if( check_if_string_is_a_number_from_1_to_100( expanded_input_filename.second[j] )) {
					add_input_and_corresponding_output_filename( expanded_input_filename.first[j], output_filename, expanded_input_filename.second[j] ) ;
				}
			}
		}
		else {
			add_input_and_corresponding_output_filename( input_filename, output_filename, "" ) ;
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
		std::size_t wildcard_pos = output_filename.find( '#' ) ;
		if( wildcard_pos != std::string::npos ) {
			expanded_output_filename.replace( wildcard_pos, 1, input_filename_wildcard_match ) ;			
		}
		m_gen_input_filenames.push_back( input_filename ) ;
		m_gen_output_filenames.push_back( expanded_output_filename ) ;
	}

	void open_gen_input_files() {
		#ifdef HAVE_BOOST_TIMER
			boost::timer timer ;
		#endif

		m_input_chain.reset( new genfile::SNPDataSourceChain() ) ;
		m_input_chain->set_moved_to_next_source_callback( boost::bind( &GenConvertProcessor::move_to_next_output_file, this, _1 )) ;
		
		std::size_t number_of_snps = 0 ;
		for( std::size_t i = 0; i < m_gen_input_filenames.size(); ++i ) {
			add_gen_file_to_chain( *m_input_chain, m_gen_input_filenames[i] ) ;
			m_gen_file_snp_counts.push_back( m_input_chain->total_number_of_snps() - number_of_snps ) ;
			number_of_snps = m_input_chain->total_number_of_snps() ;
		}
		
		m_number_of_samples_from_gen_file = m_input_chain->number_of_samples() ;
		
		#ifdef HAVE_BOOST_TIMER
			if( timer.elapsed() > 1.0 ) {
				m_cout << "Opened " << m_gen_input_filenames.size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
			}
		#endif
	}

	void add_gen_file_to_chain( genfile::SNPDataSourceChain& chain, std::string const& filename ) {
		#ifdef HAVE_BOOST_TIMER
			boost::timer file_timer ;
		#endif
			m_cout << "(Opening gen file \"" << filename << "\"...)" << std::flush ;
			try {
				std::auto_ptr< genfile::SNPDataSource > snp_data_source( genfile::SNPDataSource::create( filename )) ;
				chain.add_source( snp_data_source ) ;
			}
			catch ( genfile::FileHasTwoConsecutiveNewlinesError const& e ) {
				std::cerr << "\n!!ERROR: a GEN file was specified having two consecutive newlines.\n"
					<< "!! NOTE: popular editors, such as vim and nano, automatically add an extra newline to the file (which you can't see).\n"
					<< "!!     : Please check that each SNP in the file is terminated by a single newline.\n" ;
				throw ;
			}
		#ifdef HAVE_BOOST_TIMER
			m_cout << " (" << file_timer.elapsed() << "s)\n" ;
		#else
			m_cout << "\n" ;
		#endif			
	}

	void open_gen_output_files() {
		m_output_chain.reset( new genfile::SNPDataSinkChain() ) ;
		std::string output_filename = "" ;
		for( std::size_t i = 0 ; i < m_gen_output_filenames.size(); ++i ) {
			if( m_gen_output_filenames[i] != output_filename ) {
				output_filename = m_gen_output_filenames[i] ;
				add_gen_file_to_chain( *m_output_chain, output_filename ) ;
			}
		}
		// Set up the current output filename so we can track changes to the filename.
		
		m_current_output_filename = m_gen_output_filenames.front() ;
	}

	void add_gen_file_to_chain( genfile::SNPDataSinkChain& chain, std::string const& filename ) {
			std::auto_ptr< genfile::SNPDataSink > snp_data_sink( genfile::SNPDataSink::create( filename )) ;
			chain.add_sink( snp_data_sink ) ;
	}

	void move_to_next_output_file( std::size_t index ) {
		if( index < m_gen_input_filenames.size() ) {
			m_output_file_index = index ;
			if( m_gen_output_filenames[ index ] != m_current_output_filename ) {
				m_current_output_filename = m_gen_output_filenames[ index ] ;
				m_output_chain->move_to_next_sink() ;
			}
		}
		else {
			m_cout << "\n" ;
		}
	}

public:
	
	void write_banner( std::ostream& oStream ) const {
		oStream << "\nWelcome to " << globals::program_name << "\n"
		 	<< "(C) 2009 University of Oxford\n\n";
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
					<< "\"" << std::setw(20) << m_gen_input_filenames[i] << "\""
					<< " -> \"" << m_gen_output_filenames[i] << "\"\n";
			}
			oStream << std::string( 30, ' ' ) << "  (total " << m_input_chain->total_number_of_snps() << " snps).\n\n" ;

			if( !m_errors.empty() ) {
				for( std::size_t i = 0; i < m_errors.size(); ++i ) {
					oStream << "!! ERROR: " << m_errors[i] << "\n\n" ;
				}
				oStream << "!! Please correct the above errors and re-run " << globals::program_name << ".\n\n" ;
				throw GenConvertUsageError() ;
			}

			oStream << std::string( 72, '=' ) << "\n\n" ;
		}
		catch (...) {
			oStream << std::string( 72, '=' ) << "\n\n" ;
			throw ;
		}
	}

	void write_postamble( std::ostream& oStream ) const {
		oStream << "\n"
			<< "Thank you for using " << globals::program_name << ".\n" ;
	}

	void do_checks() {
		for( std::size_t i = 0; i < m_gen_input_filenames.size(); ++i ) {
			for( std::size_t j = 0; j < m_gen_output_filenames.size(); ++j ) {
				if( strings_are_nonempty_and_equal( m_gen_input_filenames[i], m_gen_output_filenames[j] )) {
					m_errors.push_back( "Output GEN file \"" + m_gen_output_filenames[j] +"\" also specified as input GEN file." ) ;
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
#if HAVE_BOOST_TIMER
		boost::timer timer ;
#endif
		open_gen_output_files() ;

		GenRow row ;
		row.set_number_of_samples( m_number_of_samples_from_gen_file ) ;

		m_cout << "Converting GEN files...\n" ;

		double last_time = -5.0 ;
		std::size_t number_of_snps_processed = 0 ;
		while( read_snp(row) ) {
			++number_of_snps_processed ;
			// print a message every 5 seconds.

#if HAVE_BOOST_TIMER
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
#else
			if( number_of_snps_processed % 1000 == 0 ) {
				m_cout << "Processed " << number_of_snps_processed << " SNPs (" << timer.elapsed() << "s)\n" ;
			}
#endif
			write_snp(row) ;
		}

	#if HAVE_BOOST_TIMER
		std::cerr << "Converted GEN file(s) (" << number_of_snps_processed << " SNPs) in " << timer.elapsed() << " seconds.\n" ;
	#endif
	
		m_cout << "Post-processing (updating file header, compression)..." << std::flush ;
		timer.restart() ;
		// Close the output gen file(s) now
		close_all_files() ;
	#if HAVE_BOOST_TIMER
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
		GenConvertProcessor::declare_options( options ) ;
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
		GenConvertProcessor processor( options ) ;
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
