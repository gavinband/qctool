/*
 * This program, gen-case-control-test
 */

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <numeric>
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
	std::string const program_name = "gen-compare" ;
}

struct GenCompareProcessorException: public GToolException
{
	char const* what() const throw() { return "GenCompareProcessorException" ; }
} ;

struct GenCompareUsageError: public GenCompareProcessorException
{
	char const* what() const throw() { return "GenCompareUsageError" ; }
} ;

// thrown to indicate that the numbers of files specified for -cases or -controls
// is not two
struct GenCompareFileCountError: public GenCompareProcessorException
{
	char const* what() const throw() { return "GenCompareFileCountError" ; }
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
	    options[ "-g1" ]
	        .set_description( "Path of gen file to input" )
			.set_is_required()
			.set_takes_single_value()
			.set_takes_value_by_position( 1 ) ;

	    options[ "-g2" ]
	        .set_description( "Path of gen file to input" )
			.set_is_required()
			.set_takes_single_value()
			.set_takes_value_by_position( 2 ) ;

		options [ "--force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
	}

	GenCompareProcessor( OptionProcessor const& options )
		: m_cout( std::cout.rdbuf() ),
		  m_options( options )
	{
		write_start_banner( m_cout ) ;
		setup() ;
	}

	~GenCompareProcessor() {
		write_end_banner( m_cout ) ;
	}
	
private:
	void setup() {
		try {
			get_required_filenames() ;
			open_gen_files() ;
		}
		catch( genfile::FileStructureInvalidError const& ) {
			m_cout << "\nThe GEN file specified does not seem to be in a valid gen or bgen format.\n" ;
			throw ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_cout << "\nThe GEN files specified did not all have the same sample size.\n" ;
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
		std::string case_filename, control_filename ;
		control_filename = m_options.get_value< std::string >( "-g1" ) ;
		expand_and_add_filename( &m_control_gen_filenames, control_filename ) ;
		case_filename = m_options.get_value< std::string >( "-g2" ) ;
		expand_and_add_filename( &m_case_gen_filenames, case_filename ) ;
	}

	void expand_and_add_filename( std::map< int, std::string >* filename_list_ptr, std::string const& filename ) {
		bool input_file_has_wildcard = ( filename.find( '#' ) != std::string::npos ) ;
		if( input_file_has_wildcard ) {
			std::pair< std::vector< std::string >, std::vector< std::string > >
				expanded_filename = find_files_matching_path_with_wildcard( filename, '#' ) ;

			// we only use matching filenames if the match is a number from 1 to 100
			// For such filenames, we add the filename to our list for cases.
			for( std::size_t j = 0; j < expanded_filename.first.size(); ++j ) {
				int file_number ;
				if( parse_number_from_1_to_100( expanded_filename.second[j], &file_number )) {
					add_filename( filename_list_ptr, file_number, expanded_filename.first[j] ) ;
				}
			}
		}
		else {
			add_filename( filename_list_ptr, 1, filename ) ;
		}
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

	void add_filename( std::map< int, std::string >* filename_list_ptr, int file_number, std::string const& filename ) {
		filename_list_ptr->insert( std::make_pair( file_number, filename )) ;
	}

	void open_gen_files() {
		Timer timer ;
		open_gen_files( m_control_gen_filenames, m_control_gen_input_chain ) ;
		open_gen_files( m_case_gen_filenames, m_case_gen_input_chain ) ;

		if( timer.elapsed() > 1.0 ) {
			m_cout << "Opened " << m_control_gen_filenames.size() + m_case_gen_filenames.size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
		}
	}
	
	void open_gen_files( std::map< int, std::string > const& filenames, std::auto_ptr< genfile::SNPDataSourceChain >& chain ) {
		chain.reset( new genfile::SNPDataSourceChain() ) ;
		std::map< int, std::string >::const_iterator
			i = filenames.begin(),
			end_i = filenames.end() ;
		for( ; i != end_i ; ++i ) {
			add_gen_file_to_chain( *chain, i->second ) ;
		}
	}

	void add_gen_file_to_chain( genfile::SNPDataSourceChain& chain, std::string const& filename ) {
		Timer timer ;
		m_cout << "(Opening gen file \"" << filename << "\"...)" << std::flush ;
		std::auto_ptr< genfile::SNPDataSource > snp_data_source( genfile::SNPDataSource::create( filename )) ;
		chain.add_source( snp_data_source ) ;
		m_cout << " (" << timer.elapsed() << "s)\n" ;
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
			oStream << "First set of GEN files:\n" ;
			print_gen_files( oStream, m_control_gen_filenames, m_control_gen_input_chain ) ;
			oStream << "Second set of GEN files:\n" ;
			print_gen_files( oStream, m_case_gen_filenames, m_case_gen_input_chain ) ;
			oStream << std::string( 72, '=' ) << "\n\n" ;
		}
		catch (...) {
			oStream << std::string( 72, '=' ) << "\n\n" ;
			throw ;
		}
	}

	void print_gen_files( std::ostream& oStream, std::map< int, std::string > const& filenames, std::auto_ptr< genfile::SNPDataSourceChain > const& chain ) const {
		std::map< int, std::string >::const_iterator
			i = filenames.begin(),
			end_i = filenames.end() ;
		std::size_t count = 0 ;
		for( ; i != end_i ; ++i, ++count ) {
			oStream
				<< "  (" << std::setw(6) << chain->number_of_snps_in_source( count ) << " snps)  "
				<< "\"" << std::setw(20) << i->second << "\"\n" ;
		}
		oStream << "  (total " << chain->total_number_of_snps() << " snps).\n\n" ;
	}

	void do_checks() {
	}
	
	void process() {
		unsafe_process() ;
	}

private:

	void unsafe_process() {
		Timer timer ;

		InternalStorageGenRow case_row, control_row ;
		control_row.set_number_of_samples( m_control_gen_input_chain->number_of_samples() ) ;
		case_row.set_number_of_samples( m_case_gen_input_chain->number_of_samples() ) ;

		m_cout << "Comparing SNPs in gen files...\n" ;

		double last_time = -5.0 ;
		std::size_t number_of_control_snps_matched = 0 ;
		bool these_ones_matched = true ;
		while( read_snp( *m_case_gen_input_chain, case_row )) {
			std::size_t number_of_matching_snps = 0 ;
			while( read_next_snp_with_specified_position( *m_control_gen_input_chain, control_row, case_row.SNP_position())) {
				++number_of_matching_snps ;
				++number_of_control_snps_matched ;
				// I think we'd expect all the SNP identifying data to be the same for the two snps.
				// If there are mismatches, print out some information about them.
				bool this_one_matched = compare_snps( case_row, control_row ) ;
				if( !this_one_matched ) {
					m_cout << "Rows "
						<< m_case_gen_input_chain->number_of_snps_read()
						<< " of first file, and "
						<< m_control_gen_input_chain->number_of_snps_read()
						<< " of second file differ (but match).\n" ;
				}
				these_ones_matched = these_ones_matched && this_one_matched ;
				
				// print a progress message every second.
				double time_now = timer.elapsed() ;
				if( (time_now - last_time >= 1.0) || (m_case_gen_input_chain->number_of_snps_read() == m_case_gen_input_chain->total_number_of_snps()) ) {
					print_progress( time_now, these_ones_matched ) ;
					last_time = time_now ;
					these_ones_matched = true ;
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

		std::cerr
			<< "\nProcessed all SNPs ("
			<< m_case_gen_input_chain->number_of_snps_read() << " SNPs in first file, "
			<< number_of_control_snps_matched << " matched SNPs in second file) in "
			<< std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
	
		close_all_files() ;
	}
	
	bool compare_snps( GenRow const& case_row, GenRow const& control_row ) {
		bool matched = true ;
		if( case_row != control_row ) {
			m_cout << "\nRows differ!\n" ;
			matched = false ;
		}
		if( case_row.SNPID() != control_row.SNPID() ) {
			m_cout << "\nMatching rows have differing SNPIDs " << case_row.SNPID() << " and " << control_row.SNPID() << ".\n" ;
			matched = false ;
		}
		if( case_row.RSID() != control_row.RSID() ) {
			m_cout << "\nMatching rows have differing RSIDs " << case_row.RSID() << " and " << control_row.RSID() << ".\n" ;
			matched = false ;
		}
		if( case_row.SNP_position() != control_row.SNP_position() ) {
			m_cout << "\nMatching rows have different positions " << case_row.SNP_position() << " and " << control_row.SNP_position() << ".\n" ;
			matched = false ;
		}
		if( case_row.first_allele() != control_row.first_allele() || case_row.second_allele() != control_row.second_allele() ) {
			m_cout 	<< "\nMatching rows have differing alleles "
					<< case_row.first_allele() << " " << case_row.second_allele()
					<< " and " << control_row.first_allele() << " " << control_row.second_allele() << ".\n" ;
				matched = false ;
		}
		return matched ;
	}

	void print_progress( double time_now, bool these_ones_matched ) {
		double case_progress = (static_cast< double >( m_case_gen_input_chain->number_of_snps_read() ) / m_case_gen_input_chain->total_number_of_snps()) ;
		double control_progress = (static_cast< double >( m_control_gen_input_chain->number_of_snps_read() ) / m_control_gen_input_chain->total_number_of_snps()) ;
		m_cout
			<< "\r"
			<< get_progress_bar( 30, case_progress, these_ones_matched )
			<< " (" << m_case_gen_input_chain->number_of_snps_read() << "/" << m_case_gen_input_chain->total_number_of_snps() << ")"
			<< " "
			<< get_progress_bar( 30, control_progress, these_ones_matched )
			<< " (" << m_control_gen_input_chain->number_of_snps_read() << "/" << m_control_gen_input_chain->total_number_of_snps() << ")"
			<< " (" << std::fixed << std::setprecision(1) << time_now << "s)"
			<< std::string( std::size_t(5), ' ' )
			<< std::flush ;
	}

	std::string get_progress_bar( std::size_t width, double progress, bool these_ones_matched ) {
		progress = std::min( std::max( progress, 0.0 ), 1.0 ) ;
		std::size_t visible_progress = progress * width ;

		if( visible_progress > 0 ) {
			return
				"["
				+ std::string( std::size_t( visible_progress - 1 ), '*' )
				+ std::string( std::size_t(1), these_ones_matched ? '*' : '!' )
				+ std::string( std::size_t( width - visible_progress ), ' ' )
				+ "]" ;
		}
		else {
			return "[" + std::string( std::size_t( width ), ' ' ) + "]" ;
		}
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

	bool read_next_snp_with_specified_position( genfile::SNPDataSource& snp_data_source, GenRow& row, uint32_t specified_SNP_position ) {
		return snp_data_source.read_next_snp_with_specified_position(
			boost::bind< void >( &GenRow::set_number_of_samples, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNPID, &row, _1 ),
			boost::bind< void >( &GenRow::set_RSID, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNP_position, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele1, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele2, &row, _1 ),
			boost::bind< void >( &GenRow::set_genotype_probabilities, &row, _1, _2, _3, _4 ),
			specified_SNP_position 
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
	
	std::map< int, std::string > m_case_gen_filenames ;
	std::map< int, std::string > m_control_gen_filenames ;

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
        std::cerr << "Usage: " << globals::program_name << " [options]\n"
                << options
                << "\n" ;
        return -1 ;
    }

	// main section

	try	{
		GenCompareProcessor processor( options ) ;
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
