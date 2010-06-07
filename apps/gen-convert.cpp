/*
 * This program, gen-convert, compares two gen files (or numbered collections of gen files)
 * on a SNP-by-SNP basis.
 */

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <numeric>
#include <boost/bind.hpp>
#include "Timer.hpp"
#include "GenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "OptionProcessor.hpp"
#include "FileUtil.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/SNPDataSinkChain.hpp"
#include "string_utils/string_utils.hpp"
#include "string_utils/parse_utils.hpp"
#include "wildcard.hpp"
#include "InputToOutputFilenameMapper.hpp"

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

// thrown to indicate that a specified GEN file could not be found.
struct GenConvertNoGenFileMatchesFound: public GenConvertProcessorException
{
	GenConvertNoGenFileMatchesFound( std::string const& filename ): m_filename( filename ) {}
	~GenConvertNoGenFileMatchesFound() throw() {}

	char const* what() const throw() { return "GenConvertNoGenFileMatchesFound" ; }
	
	std::string const& filename() const { return m_filename ; }
private:
	std::string const m_filename ;
} ;

struct ChromosomeNotSuppliedError: public GenConvertProcessorException
{
	char const* what() const throw() { return "ChromosomeNotSuppliedError" ; }
} ;

struct GenConvertProcessor
{
public:
	static void declare_options( OptionProcessor & options ) {
		
		// File options
		options.declare_group( "File options" ) ;
	    options[ "-g" ]
	        .set_description( "Path of gen file to input" )
	        .set_is_required()
			.set_takes_single_value()
			.set_takes_value_by_position(1) ;

	    options[ "-og" ]
	        .set_description( "Path of gen file to output" )
			.set_takes_single_value()
			.set_takes_value_by_position(2) ;

		options.declare_group( "Other options" ) ;
			
		options [ "-force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
			
		options [ "-chromosome" ]
			.set_description( "Override the chromosomal values inferred from the input file(s).")
			.set_takes_single_value() ;
		options [ "-sort" ]
			.set_description( "Sort the output files by SNP chromosome, position, and RSID.  This is only "
				"supported if bgen output format is used." ) ;
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
			get_chromosome_if_specified() ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_cout << "The GEN files specified did not all have the same sample size.\n" ;
			throw ;
		}
		catch( GenConvertFileCountMismatchError const& ) {
			m_cout << "The number of files specified for -g must match that for -og.\n" ;
			throw ;
		}
		catch( GenConvertFileWildcardMismatchError const& ) {
			m_cout << "The output file contained a wildcard character, so the input should also.\n" ;
			throw ;
		}
		catch( ChromosomeNotSuppliedError const& ) {
			m_cout << string_utils::wrap( "Error: the input filename contained no wildcard (#) character.  The gen file format does "
				"not support chromosome entries, and I can't deduce it from context, so you must specify the "
				"chromosome using the -chromosome option.\n", 80 );
			throw ;
		}
	}

	void get_required_filenames() {
		std::vector< std::string > gen_filenames, gen_output_filenames ;
		if( m_options.check_if_option_was_supplied( "-g" ) ) {
			gen_filenames = m_options.get_values< std::string >( "-g" ) ;
		}
		if( m_options.check_if_option_was_supplied( "-og" ) ) {
			gen_output_filenames = m_options.get_values< std::string >( "-og" ) ;
		} else {
			gen_output_filenames.resize( gen_filenames.size() ) ;
		}

		if( gen_filenames.size() != gen_output_filenames.size() ) {
			throw GenConvertFileCountMismatchError() ;
		}

		for( std::size_t i = 0; i < gen_filenames.size(); ++i ) {
			if(( gen_filenames[i].find( '#' ) == std::string::npos )
				&& genfile::filename_indicates_gen_format( gen_filenames[i] )
				&& ( !m_options.check_if_option_was_supplied( "-chromosome" ))) {
					throw ChromosomeNotSuppliedError() ;
			}
		}

		m_gen_file_mapper.add_filename_pairs( gen_filenames, gen_output_filenames ) ;
	}

	void open_gen_input_files() {
		Timer timer ;

		m_input_chain.reset( new genfile::SNPDataSourceChain() ) ;
		m_input_chain->set_moved_to_next_source_callback( boost::bind( &GenConvertProcessor::move_to_next_output_file, this, _1 )) ;
		
		for( std::size_t i = 0; i < m_gen_file_mapper.input_files().size(); ++i ) {
			add_gen_file_to_chain( *m_input_chain, m_gen_file_mapper.input_file(i), m_gen_file_mapper.matched_part(i) ) ;
		}
		
		m_number_of_samples_from_gen_file = m_input_chain->number_of_samples() ;
		
		if( timer.elapsed() > 1.0 ) {
			m_cout << "Opened " << m_gen_file_mapper.input_files().size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
		}
	}

	void add_gen_file_to_chain( genfile::SNPDataSourceChain& chain, std::string const& filename, std::string const& matched_part ) {
		boost::timer file_timer ;
		m_cout << "(Opening gen file \"" << filename << "\"...)" << std::flush ;
		try {
			std::auto_ptr< genfile::SNPDataSource > snp_data_source( genfile::SNPDataSource::create( filename, matched_part )) ;
			chain.add_source( snp_data_source ) ;
		}
		catch ( genfile::FileHasTwoConsecutiveNewlinesError const& e ) {
			std::cerr << "\n!!ERROR: a GEN file was specified having two consecutive newlines.\n"
				<< "!! NOTE: popular editors, such as vim and nano, automatically add an extra newline to the file (which you can't see).\n"
				<< "!!     : Please check that each SNP in the file is terminated by a single newline.\n" ;
			throw ;
		}
		m_cout << " (" << file_timer.elapsed() << "s)\n" ;
	}

	void open_gen_output_files() {
		m_output_chain.reset( new genfile::SNPDataSinkChain() ) ;
		std::string output_filename = "" ;
		for( std::size_t i = 0 ; i < m_gen_file_mapper.output_filenames().size(); ++i ) {
				add_gen_file_to_chain( *m_output_chain, m_gen_file_mapper.output_filename(i) ) ;
		}

		// Set up the current output filename so we can track changes to the filename.
		if( m_gen_file_mapper.output_filenames().size() > 0 ) {
			m_current_output_filename = m_gen_file_mapper.output_filenames().front() ;
		}
	}

	void add_gen_file_to_chain( genfile::SNPDataSinkChain& chain, std::string const& filename ) {
			std::auto_ptr< genfile::SNPDataSink > snp_data_sink(
				genfile::SNPDataSink::create(
					filename,
					"Created by gen-convert",
					m_options.check_if_option_was_supplied( "-sort" )
				)
			) ;
			chain.add_sink( snp_data_sink ) ;
	}

	void move_to_next_output_file( std::size_t index ) {
		if( index < m_gen_file_mapper.input_files().size() ) {
			if( m_gen_file_mapper.output_filenames().size() > 0 ) {
				if( m_gen_file_mapper.filename_corresponding_to( index ) != m_output_chain->index_of_current_sink() ) {
					m_output_chain->move_to_next_sink() ;
				}
			}
		}
	}
	
	void get_chromosome_if_specified() {
		if( m_options.check_if_option_was_supplied( "-chromosome" )) {
			m_override_chromosome = true ;
			m_chromosome = m_options.get_value< std::string >( "-chromosome" ) ;
		}
		else {
			m_override_chromosome = false ;
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
			oStream << std::setw(30) << "Input and output GEN files:\n" ;
			for( std::size_t i = 0; i < m_gen_file_mapper.input_files().size(); ++i ) {
				oStream
					<< "  - "
					<< " (" << std::setw(6) << (m_input_chain->number_of_snps_in_source(i)) << " snps)  "
					<< "\"" << std::setw(20) << m_gen_file_mapper.input_file(i) << "\""
					<< " -> \"" << m_gen_file_mapper.output_filename( m_gen_file_mapper.filename_corresponding_to( i ) ) << "\"\n";
			}
			oStream << "  - (total " << m_input_chain->total_number_of_snps() << " snps).\n\n" ;

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
		for( std::size_t i = 0; i < m_gen_file_mapper.input_files().size(); ++i ) {
			for( std::size_t j = 0; j < m_gen_file_mapper.output_filenames().size(); ++j ) {
				if( strings_are_nonempty_and_equal( m_gen_file_mapper.input_file(i), m_gen_file_mapper.output_filenames()[j] )) {
					m_errors.push_back( "Output GEN file \"" + m_gen_file_mapper.output_filenames()[j] +"\" also specified as input GEN file." ) ;
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
		Timer timer ;
		open_gen_output_files() ;

		InternalStorageGenRow row ;
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

			if( m_override_chromosome ) {
				row.set_chromosome( m_chromosome ) ;
			}
			
			write_snp(row) ;
		}

		m_cout << "\nConverted GEN file(s) (" << number_of_snps_processed << " SNPs) in " << timer.elapsed() << " seconds.\n" ;
	
		timer.restart() ;
		// Close the output gen file(s) now
		close_all_files() ;
		std::cerr << " (" << timer.elapsed() << "s)\n" ;
	}

	bool read_snp( GenRow& row ) {
		return row.read_from_source( *m_input_chain ) ;
	}

	bool write_snp( GenRow& row ) {
		return row.write_to_sink( *m_output_chain ) ;
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
	
	std::size_t m_number_of_samples_from_gen_file ;
	
	InputToOutputFilenameMapper m_gen_file_mapper ;

	std::size_t m_output_file_index ;

	std::vector< std::string > m_errors ;
	
	bool m_override_chromosome ;
	genfile::Chromosome m_chromosome ;
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
