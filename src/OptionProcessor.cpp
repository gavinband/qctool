#include <string>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "GenRow.hpp"
#include "GToolException.hpp"
#include "OptionDefinition.hpp"
#include "OptionProcessor.hpp"

OptionProcessingException::OptionProcessingException( std::string option, std::vector< std::string > values, std::string msg )
: 	m_option( option ),
	m_values( values ),
	m_msg( msg )
{}

OptionProcessingException::OptionProcessingException( std::string option, std::string msg )
: 	m_option( option ),
	m_msg( msg )
{}

OptionProcessingException::~OptionProcessingException() throw()
{} ;

OptionParseException::OptionParseException( std::string option, std::vector< std::string > values, std::string msg )
: OptionProcessingException( option, values, msg )
{}

OptionValueInvalidException::OptionValueInvalidException( std::string option, std::vector< std::string > values, std::string msg )
: OptionProcessingException( option, values, msg )
{}


char const* OptionProcessingException::OptionProcessingException::what() const throw() {
	try {
		std::ostringstream ostr ;
		if( m_option.size() > 0 ) {
			ostr << "For option \"" << m_option << "\": " ;
		}
		ostr
			<< m_msg
			<< " (values supplied:" ;
		for( std::size_t i = 0; i < m_values.size() ; ++i ) {
			if( i > 0 )
				ostr << "," ;
			ostr << m_values[i] ;
		}
		ostr << ")" ;
		set_message( ostr.str() ) ;
		return GToolException::what() ;
	}
	catch (... ) {
		return "OptionProcessingException::what(): exception formatting string." ;
	}
}


std::ostream& operator<<( std::ostream& aStream, OptionProcessor::OptionDefinitions const& option_definitions ) {
    // Find longest option and longest option description
    OptionProcessor::OptionDefinitions::const_iterator
        defn_i = option_definitions.begin(),
        defn_end = option_definitions.end() ;

    std::size_t max_option_length = 0, max_description_length = 0 ;
     
    for( ; defn_i != defn_end; ++defn_i ) {
        max_option_length = std::max( max_option_length, defn_i->first.size()) ;
        max_description_length = std::max( max_description_length, defn_i->second.description().size() ) ;
    }

	// print options
    defn_i = option_definitions.begin() ;
    for( ; defn_i != defn_end; ++defn_i ) {
        aStream << std::setw( max_option_length+1 )
                << std::right
                << defn_i->first
				<< "  : "
				<< defn_i->second.description() ;
		
		aStream << ".\n" ;
    }

    return aStream ;
}


OptionProcessor::OptionProcessor() {}
OptionProcessor::OptionProcessor( OptionProcessor const& other )
	: m_option_definitions( other.m_option_definitions ),
		m_option_values( other.m_option_values ),
		m_unknown_options( other.m_unknown_options )
{}

OptionProcessor::~OptionProcessor() {} ;
		
OptionProcessor::OptionDefinitions const& OptionProcessor::option_definitions() const { return m_option_definitions ; }

OptionDefinition& OptionProcessor::operator[]( std::string const& arg ) {
	return m_option_definitions[ arg ] ;
}

OptionDefinition const& OptionProcessor::operator[] ( std::string const& arg ) const {
	OptionDefinitions::const_iterator defn_i = m_option_definitions.find( arg ) ;
	if( defn_i == m_option_definitions.end() )
		throw OptionValueInvalidException( arg, std::vector< std::string >(), "Option \"" + arg + "\" is not defined." ) ;
	return defn_i->second ;
}

// Parse the options from argv, performing all needed checks.
void OptionProcessor::process( int argc, char** argv ) {
	parse_options( argc, argv ) ;
	check_required_options_are_supplied() ;
	preprocess_option_values() ;
	check_option_values() ;
}

// Parse the options.  Store option values.  Ignore, but store unknown args for later reference
void OptionProcessor::parse_options( int argc, char** argv ) {
	int i = 0;
	while( ++i < argc ) {
		std::string arg = argv[i] ;
		std::map< std::string, OptionDefinition >::const_iterator arg_i = m_option_definitions.find( arg ) ;
		if( !try_to_parse_named_option_and_values( argc, argv, i )) {
			if( !try_to_parse_positional_option( argc, argv, i )) {
				m_unknown_options[i] = argv[i] ;
			} ;
		}
	}
	
	if( !m_unknown_options.empty()) {
		process_unknown_options() ;
	}
}	

bool OptionProcessor::try_to_parse_named_option_and_values( int argc, char** argv, int i ) {
	std::string arg = argv[i] ;
	std::map< std::string, OptionDefinition >::const_iterator arg_i = m_option_definitions.find( arg ) ;
	if( arg_i != m_option_definitions.end() ) {
		if( arg_i->second.number_of_values_per_use() > 0 ) {
			for( unsigned int value_i = 0; value_i < arg_i->second.number_of_values_per_use(); ++value_i ) {
				if( ++i == argc ) {
					std::ostringstream ostr ;
					ostr << "Option\"" << arg << "\" takes " << arg_i->second.number_of_values_per_use() << " values per use." ;
					throw OptionParseException( arg, m_option_values[ arg ], ostr.str() ) ;
				}
				else {
					m_option_values[ arg ].push_back( argv[i] ) ;
				}
			}
		}
		else {
			m_option_values[ arg ].push_back( "1" ) ;
		}
		return true ;
	}
	return false ;
}

bool OptionProcessor::try_to_parse_positional_option( int argc, char** argv, int position ) {
	OptionDefinitions::const_iterator
		defn_i = m_option_definitions.begin(),
		end_defn_i = m_option_definitions.end() ;
		
	for( ; defn_i != end_defn_i; ++defn_i ) {
		if( defn_i->second.takes_value_by_position() && defn_i->second.position() == position ) {
			m_option_values[ defn_i->first ].push_back( argv[ position ] ) ;
			return true ;
		}
	}
	return false ;
}

void OptionProcessor::process_unknown_options() {
	std::string unknown_args ;
	std::map< int, std::string >::const_iterator i = m_unknown_options.begin() ;
	for( ; i != m_unknown_options.end(); ++i ) {
		if( i != m_unknown_options.begin() )
			unknown_args += ", " ;
		unknown_args += i->second ;
	}
	throw OptionParseException( "", std::vector< std::string >(), "The following options were not recognised: " + unknown_args + "." ) ;
}


void OptionProcessor::check_required_options_are_supplied() {
	std::map< std::string, OptionDefinition >::const_iterator
		defn_i = m_option_definitions.begin(), 
		defn_end = m_option_definitions.end() ;

	for( ; defn_i != defn_end; ++defn_i ) {
		if( defn_i->second.is_required() && !check_if_option_was_supplied( defn_i->first ) && !defn_i->second.has_default_value()) {
			throw OptionProcessingException( defn_i->first, std::vector< std::string >(), "Option \"" + defn_i->first + "\" must be supplied." ) ;
		}
	}
}

void OptionProcessor::preprocess_option_values() {
	OptionValues::iterator
		arg_i = m_option_values.begin(),
		arg_end = m_option_values.end() ;

	for( ; arg_i != arg_end; ++arg_i ) {
		OptionDefinitions::const_iterator defn_i
			= m_option_definitions.find( arg_i->first ) ;
	
		arg_i->second = defn_i->second.preprocess_option_values( arg_i->first, arg_i->second ) ;
	}
}

void OptionProcessor::check_option_values() {
	OptionValues::const_iterator
		arg_i = m_option_values.begin(),
		arg_end = m_option_values.end() ;

	for( ; arg_i != arg_end; ++arg_i ) {
		OptionDefinitions::const_iterator defn_i
			= m_option_definitions.find( arg_i->first ) ;

		// Run user-supplied checker if supplied
	   defn_i->second.check_option_values( arg_i->first, arg_i->second ) ; 
	}
}

// check if the given option (which must be valid) was supplied.
bool OptionProcessor::check_if_option_was_supplied( std::string const& arg ) const {
	return ( m_option_values.find( arg ) != m_option_values.end() ) ;
}

std::string OptionProcessor::get_default_value( std::string const& arg ) const {
	std::map< std::string, OptionDefinition >::const_iterator defn_i
		= m_option_definitions.find( arg ) ;
	assert( defn_i != m_option_definitions.end() ) ;
	assert( defn_i->second.has_default_value() ) ;
	return defn_i->second.default_value() ;
}

std::ostream& operator<<( std::ostream& aStream, OptionProcessor const& options ) {
	return aStream << "Options:\n" << options.option_definitions() ;
}

// get the values of the given option as a whitespace-separated string
// Note: from this function, it may not possible to determine if an option value
// contained whitespace or if several values were supplied.  Use the std::vector form
// below to avoid ambiguity.
template<>
std::string OptionProcessor::get_value( std::string const& arg ) const {
	OptionValues::const_iterator arg_i ;
	arg_i = m_option_values.find( arg ) ;
	
	if( arg_i == m_option_values.end() ) {
		return get_default_value( arg ) ;
	}

	assert( arg_i->second.size() == 1 ) ;
	return arg_i->second[0] ;
}

// get the value(s) of the given option as a vector of strings.
template<>
std::vector< std::string > OptionProcessor::get_values< std::string >( std::string const& arg ) const {
	OptionValues::const_iterator arg_i ;
	arg_i = m_option_values.find( arg ) ;
	
	if( arg_i == m_option_values.end() ) {
		return std::vector< std::string > ( 1, get_default_value( arg )) ;
	}

	return arg_i->second ;
}


