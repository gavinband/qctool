
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include <set>
#include <iostream>
#include <boost/format.hpp>
#include <boost/optional.hpp>

#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/ProgramFlow.hpp"

#include "genfile/Error.hpp"
#include "genfile/GenomePositionRange.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/FileUtils.hpp"

namespace globals {
	std::string const program_name = "binit" ;
	std::string const program_version = "0.1" ;
}

struct BinitOptions: public appcontext::CmdLineOptionProcessor {
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;

		options[ "-bin-size" ]
			.set_description( "Number of input items to put in each bin. "
				"A value of 0 specifies a single bin containing all the data." )
			.set_takes_single_value()
			.set_is_required()
			.set_default_value( 0 )
		;

		options[ "-assume-range" ]
			.set_description( "A chromosome and position range (e.g. <chr>:<start>-<end>). "
				"If this is given it is assumed that the N numbers passed into the standard input "
				"correspond to the positions start - end on the given chromosome (and binit "
				"will check that N = end - start + 1)."
			)
			.set_takes_single_value()
		;

		options[ "-at-positions" ]
			.set_description( "Specify a file containing a list of positions. "
				"Only values at the specified positions will be used. "
				"If -assume-range is specified, positions should be in the form chr:position. "
				"Otherwise, they should be indices in the range 1...N"
			)
			.set_takes_single_value()
		;

		options[ "-ignore-lines-starting-with" ]
			.set_description( "Specify a string to treat as a comment, e.g. such that lines beginning "
				"with this string will be ignored.  Multiple values may be given." )
			.set_takes_values_until_next_option()
			.set_maximum_multiplicity( 100 )
		;

		options[ "-analysis-name" ]
			.set_description( "Specify a name for the analysis. This will be included in the output." )
			.set_takes_single_value()
			.set_default_value( appcontext::get_current_time_as_string() ) ;
		;

		options[ "-precision" ]
			.set_description( "Number of significant digits to output results to." )
			.set_takes_single_value()
			.set_default_value( 5 )
		;

		options[ "-log" ]
			.set_description( "Path of log file." )
			.set_takes_single_value()
		;
	}
} ;

struct BinitApplication: public appcontext::ApplicationContext {
public:
	BinitApplication( int argc, char **argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version,
			std::auto_ptr< appcontext::OptionProcessor >( new BinitOptions ),
			argc,
			argv,
			"-log"
		),
		m_bin_size( options().get< std::size_t >( "-bin-size" ) )
	{
		if( options().check( "-assume-range" )) {
			m_range = genfile::GenomePositionRange::parse( options().get< std::string >( "-assume-range" ) ) ;
		}
		if( options().check( "-ignore-lines-starting-with" )) {
			m_ignore_strings = options().get_values< std::string >( "-ignore-lines-starting-with" ) ;
		}
		if( options().check( "-at-positions" )) {
			m_inclusions = std::set< genfile::VariantEntry >() ;
			std::auto_ptr< std::istream > stream( genfile::open_text_file_for_input( options().get< std::string >( "-at-positions" )) ) ;
			std::string line ;
			while( std::getline( *stream, line ) ) {
				if( m_range ) {
					m_inclusions->insert( genfile::GenomePosition( line )) ;
				} else {
					m_inclusions->insert( genfile::string_utils::to_repr< int64_t >( line )) ;
				}
			}
		}
	}

	void run() {
		process( std::cin ) ;
	}
private:
	
	std::size_t const m_bin_size ;
	std::vector< std::string > m_ignore_strings ;
	boost::optional< genfile::GenomePositionRange > m_range ;
	boost::optional< std::set< genfile::VariantEntry > > m_inclusions ;
	
private:
	void process( std::istream& input ) {
		try {
			unsafe_process( input ) ;
		} catch( genfile::InputError const& e ) {
			std::cerr << "!! (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		} catch( genfile::string_utils::StringConversionError const& e ) {
			std::cerr << "!! (" << e.what() << "): input should be numerical.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}
	
	void unsafe_process( std::istream& input ) {
		boost::optional< std::string > analysis_name ;
		if( options().check( "-analysis-name" ) ) {
			analysis_name = options().get< std::string >( "-analysis-name" ) ;
			std::cout << "analysis\t" ;
		}
		if( m_range ) {
			std::cout << "chromosome\tstart\tend\tcount\taverage\n" ;
		} else {
			std::cout << "start\tend\tcount\taverage\n" ;
		}

		char const tab = '\t' ;
		std::cout << std::setprecision( options().get< std::size_t >( "-precision" ) ) ;

		std::size_t count = 0 ;
		std::size_t index = 0 ;
		double accumulation = 0 ;
		std::string line ;
		while( input >> line ) {
			bool ignore = false ;
			if( m_ignore_strings.size() > 0 ) {
				for( std::size_t i = 0; !ignore && i < m_ignore_strings.size(); ++i ) {
					if( line.compare( 0, m_ignore_strings[i].size(), m_ignore_strings[i] ) == 0 ) {
						// ignore a line.
						std::getline( input, line ) ;
						ignore = true ;
					}
				}
			}

			if( !ignore ) {
				++index ;
			}

			if( m_inclusions && !ignore ) {
				// only include if it's in the inclusions file
				if( m_range ) {
					genfile::GenomePosition position( m_range->chromosome(), m_range->start().position() + index ) ;
					ignore = ( m_inclusions->find( position ) == m_inclusions->end() ) ;
				} else {
					int64_t position = index ;
					ignore = ( m_inclusions->find( position ) == m_inclusions->end() ) ;
				}
			}
			
			if( !ignore ) {
				accumulation += genfile::string_utils::to_repr< double >( line ) ;
				++count ;
			}
			
			if( m_bin_size > 0 && ( index % m_bin_size ) == 0 ) {
				genfile::VariantEntry mean ;
				if( count > 0 ) {
					mean = ( accumulation / count ) ;
				}
				if( analysis_name ) {
					std::cout << (*analysis_name) << tab ;
				}
				if( m_range ) {
					std::cout
						<< m_range->chromosome()
						<< tab
						<< (m_range->start().position() + index - m_bin_size )
						<< tab
						<< (m_range->start().position() + index - 1 )
						<< tab
						<< count
						<< tab
						<< mean
						<< "\n" ;
				} else {
					std::cout
						<< index - m_bin_size + 1
						<< tab
						<< index
						<< tab
						<< count
						<< tab
						<< mean
						<< "\n" ;
				}
				accumulation = 0 ;
				count = 0 ;
			}
		}
		std::size_t const last_bin_count = ( m_bin_size > 0 ) ? (count % m_bin_size) : count ;
		if( last_bin_count != 0 ) {
			double const mean = ( accumulation / last_bin_count ) ;
			if( analysis_name ) {
				std::cout << (*analysis_name) << tab ;
			}
			if( m_range ) {
				std::cout
					<< m_range->chromosome()
					<< tab
					<< (m_range->start().position() + index - ( index % m_bin_size ) )
					<< tab
					<< (m_range->start().position() + index - 1 )
					<< tab
					<< last_bin_count
					<< tab
					<< mean
					<< "\n" ;
			} else {
				std::cout
					<< index - (index % m_bin_size) + 1
					<< tab
					<< index
					<< tab
					<< last_bin_count
					<< tab
					<< mean
					<< "\n" ;
			}
		}
		
		if( m_range ) {
			std::size_t expected = m_range->end().position() - m_range->start().position() + 1 ;
			if( index != expected ) {
				throw genfile::MalformedInputError(
					"(stdin)",
					( boost::format( "Length of range (%d) does not match number of elements in input (%d)." ) % expected % index ).str(),
					index
				) ;
			}
		}
	}
} ;


int main( int argc, char **argv ) {
	try {
		BinitApplication app( argc, argv ) ;	
		app.run() ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
