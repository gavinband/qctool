
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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

namespace globals {
	std::string const program_name = "binit" ;
	std::string const program_version = "0.1" ;
}

struct BinitOptions: public appcontext::CmdLineOptionProcessor {
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;

		options[ "-bin-size" ]
			.set_description( "Number of input items to put in each bin." )
			.set_takes_single_value()
			.set_is_required()
		;

		options[ "-assume-range" ]
			.set_description( "A chromosome and position range (e.g. <chr>:<start>-<end>). "
				"If this is given it is assumed that the N numbers passed into the standard input "
				"correspond to the positions start - end on the given chromosome (and binit "
				"will check that N = end - start + 1)."
			)
			.set_takes_single_value()
		;

		options[ "-ignore-lines-starting-with" ]
			.set_description( "Specify a string to treat as a comment, e.g. such that lines beginning "
				"with this string will be ignored.  Multiple values may be given." )
			.set_takes_values_until_next_option()
			.set_maximum_multiplicity( 100 )
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
		if( m_bin_size == 0 ) {
			throw genfile::BadArgumentError(
				"BinitApplication::BinitApplication()",
				"-bin-size",
				"Expected a positive number"
			) ;
		}
		if( options().check( "-assume-range" )) {
			m_range = genfile::GenomePositionRange::parse( options().get< std::string >( "-assume-range" ) ) ;
		}
		if( options().check( "-ignore-lines-starting-with" )) {
			m_ignore_strings = options().get_values< std::string >( "-ignore-lines-starting-with" ) ;
		}
	}

	void run() {
		process( std::cin ) ;
	}
private:
	
	std::size_t const m_bin_size ;
	std::vector< std::string > m_ignore_strings ;
	boost::optional< genfile::GenomePositionRange > m_range ;
	
private:
	void process( std::istream& input ) {
		try {
			unsafe_process( input ) ;
		} catch( genfile::InputError const& e ) {
			std::cerr << "!! (" << e.what() << ": " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}
	
	void unsafe_process( std::istream& input ) {
		if( m_range ) {
			std::cout << "chromosome start end count average\n" ;
		} else {
			std::cout << "start end count average\n" ;
		}
		
		std::cout << std::setprecision( options().get< std::size_t >( "-precision" ) ) ;

		std::size_t count = 0 ;
		double accumulation = 0 ;
		std::string line ;
		while( input >> line ) {
			bool ignore = false ;
			if( m_ignore_strings.size() > 0 ) {
				for( std::size_t i = 0; !ignore && i < m_ignore_strings.size(); ++i ) {
					if( line.compare( 0, m_ignore_strings[i].size(), m_ignore_strings[i] ) == 0 ) {
						ignore = true ;
					}
				}
			}
			if( !ignore ) {
				accumulation += genfile::string_utils::to_repr< double >( line ) ;
				if( ((++count) % m_bin_size) == 0 ) {
					double const mean = ( accumulation / m_bin_size ) ;
					if( m_range ) {
						std::cout
							<< m_range->chromosome()
							<< " "
							<< (m_range->start().position() + count - m_bin_size )
							<< " "
							<< (m_range->start().position() + count - 1)
							<< " "
							<< m_bin_size
							<< " "
							<< mean
							<< "\n" ;
					} else {
						std::cout
							<< count - m_bin_size + 1
							<< " "
							<< count
							<< " "
							<< m_bin_size
							<< " "
							<< mean
							<< "\n" ;
					}
					accumulation = 0 ;
				}
			}
		}
		std::size_t const last_bin_count = (count % m_bin_size) ;
		if( last_bin_count != 0 ) {
			double const mean = ( accumulation / last_bin_count ) ;
			if( m_range ) {
				std::cout
					<< m_range->chromosome()
					<< " "
					<< (m_range->start().position() + count - last_bin_count )
					<< " "
					<< (m_range->start().position() + count - 1)
					<< " "
					<< last_bin_count
					<< " "
					<< mean
					<< "\n" ;
			} else {
				std::cout
					<< (count - (count % m_bin_size) + 1 )
					<< " "
					<< count
					<< " "
					<< last_bin_count
					<< " "
					<< mean
					<< "\n" ;
			}
		}
		
		if( m_range ) {
			std::size_t expected = m_range->end().position() - m_range->start().position() + 1 ;
			if( count != expected ) {
				throw genfile::MalformedInputError(
					"(stdin)",
					( boost::format( "Length of range (%d) does not match number of elements in input (%d)." ) % expected % count ).str(),
					count
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
