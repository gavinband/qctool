#include <string>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "GToolException.hpp"
#include "OptionDefinition.hpp"
#include "OptionProcessor.hpp"
#include "string_utils.hpp"

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
		set_message( m_msg ) ;
		return GToolException::what() ;
	}
	catch (... ) {
		return "OptionProcessingException::what(): exception formatting string." ;
	}
}


OptionProcessor::OptionProcessor() {}

OptionProcessor::OptionProcessor( OptionProcessor const& other )
	: m_option_definitions( other.m_option_definitions ),
		m_option_values( other.m_option_values ),
		m_unknown_options( other.m_unknown_options ),
		m_current_group( other.m_current_group ),
		m_option_groups( other.m_option_groups ),
		m_help_option_name( other.m_help_option_name )
{}

OptionProcessor::~OptionProcessor() {} ;

OptionProcessor::OptionDefinitions const& OptionProcessor::option_definitions() const { return m_option_definitions ; }

OptionDefinition& OptionProcessor::operator[]( std::string const& arg ) {
	// Prevent modification of an option that's already defined.
	assert( m_option_definitions.find( arg ) == m_option_definitions.end() ) ;
	m_option_definitions[arg].set_group( m_current_group ) ;
	m_option_groups[ m_current_group ].insert( arg ) ;
	if( std::count( m_option_group_names.begin(), m_option_group_names.end(), m_current_group ) == 0 ) {
		m_option_group_names.push_back( m_current_group ) ;
	}
	return m_option_definitions[ arg ] ;
}

OptionDefinition const& OptionProcessor::operator[] ( std::string const& arg ) const {
	OptionDefinitions::const_iterator defn_i = m_option_definitions.find( arg ) ;
	assert( defn_i != m_option_definitions.end() ) ;
	return defn_i->second ;
}

void OptionProcessor::declare_group( std::string const& group ) {
	m_option_groups[ group ] ;
	m_current_group = group ;
}

void OptionProcessor::option_excludes_option( std::string const& excluding_option, std::string const& excluded_option ) {
	assert( m_option_definitions.find( excluding_option ) != m_option_definitions.end() ) ;
	assert( m_option_definitions.find( excluded_option ) != m_option_definitions.end() ) ;

	m_option_exclusions[excluding_option].insert( excluded_option ) ;
}


void OptionProcessor::option_excludes_group( std::string const& excluding_option, std::string const& excluded_option_group ) {
	assert( m_option_definitions.find( excluding_option ) != m_option_definitions.end() ) ;
	std::map< std::string, std::set< std::string > >::const_iterator
		option_group_i = m_option_groups.find( excluded_option_group ) ;
	assert( option_group_i != m_option_groups.end() ) ;
	std::set< std::string >::const_iterator
		i = option_group_i->second.begin(),
		end = option_group_i->second.end() ;

	for( ; i != end; ++i ) {
		option_excludes_option( excluding_option, *i ) ;
	}
}

// Parse the options from argv, performing all needed checks.
void OptionProcessor::process( int argc, char** argv ) {
	// calculate_option_groups() ;
	parse_options( argc, argv ) ;
	check_required_options_are_supplied() ;
	check_mutually_exclusive_options_are_not_supplied() ;
	check_option_values() ;
	preprocess_option_values() ;
}

void OptionProcessor::calculate_option_groups() {
	std::map< std::string, OptionDefinition >::const_iterator
		defn_i = m_option_definitions.begin(), 
		defn_end = m_option_definitions.end() ;
	for( ; defn_i != defn_end; ++defn_i ) {
		m_option_groups[ defn_i->second.group() ].insert( defn_i->first ) ;
	}
}

// Parse the options.  Store option values.  Ignore, but store unknown args for later reference
void OptionProcessor::parse_options( int argc, char** argv ) {
	int i = 0;
	while( ++i < argc ) {
		std::string arg = argv[i] ;
		if( arg == m_help_option_name ) {
			throw OptionProcessorHelpRequestedException() ;
		}
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

bool OptionProcessor::try_to_parse_named_option_and_values( int argc, char** argv, int& i ) {
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

bool OptionProcessor::try_to_parse_positional_option( int argc, char** argv, int& position ) {
	OptionDefinitions::const_iterator
		defn_i = m_option_definitions.begin(),
		end_defn_i = m_option_definitions.end() ;
		
	for( ; defn_i != end_defn_i; ++defn_i ) {
		if( defn_i->second.takes_value_by_position() && defn_i->second.position() == position ) {
			m_option_values[ defn_i->first ].push_back( argv[ position ] ) ;
			// ++position ;
			return 1 ;
		}
	}
	
	return 0 ;
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


void OptionProcessor::check_required_options_are_supplied() const {
	std::map< std::string, OptionDefinition >::const_iterator
		defn_i = m_option_definitions.begin(), 
		defn_end = m_option_definitions.end() ;

	for( ; defn_i != defn_end; ++defn_i ) {
		if( defn_i->second.is_required() && !check_if_option_was_supplied( defn_i->first ) && !defn_i->second.has_default_value()) {
			throw OptionProcessingException( defn_i->first, std::vector< std::string >(), "Option \"" + defn_i->first + "\" must be supplied." ) ;
		}
	}
}

void OptionProcessor::check_mutually_exclusive_options_are_not_supplied() const {
	std::map< std::string, std::set< std::string > > ::const_iterator
		i = m_option_exclusions.begin(),
		end_i = m_option_exclusions.end() ;

	for( ; i != end_i ; ++i ) {
		if( check_if_option_was_supplied( i->first )) {
			std::set< std::string >::const_iterator
				j = i->second.begin(),
				end_j = i->second.end() ;
			for( ; j != end_j ; ++j ) {
				if( check_if_option_was_supplied( *j )) {
					throw OptionProcessorMutuallyExclusiveOptionsSuppliedException( i->first, *j ) ;
				}
			}
		}
	}
}


void OptionProcessor::preprocess_option_values() {
	// Do nothing in base class.
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

bool OptionProcessor::check_if_option_was_supplied_in_group( std::string const& group_name ) const {
	std::map< std::string, std::set< std::string > >::const_iterator
		group_i = m_option_groups.find( group_name ) ;
	assert( group_i != m_option_groups.end() ) ;
	std::set< std::string >::const_iterator
		option_i = group_i->second.begin(),
		end_option_i = group_i->second.end() ;

	for( ; option_i != end_option_i; ++option_i ) {
		if( check_if_option_was_supplied( *option_i )) {
			return true ;
		}
	}
	return false ;
}

std::string OptionProcessor::get_default_value( std::string const& arg ) const {
	std::map< std::string, OptionDefinition >::const_iterator defn_i
		= m_option_definitions.find( arg ) ;
	assert( defn_i != m_option_definitions.end() ) ;
	assert( defn_i->second.has_default_value() ) ;
	return defn_i->second.default_value() ;
}

std::ostream& operator<<( std::ostream& aStream, OptionProcessor const& options ) {
	options.format_options( aStream ) ;
	return aStream ;
}


void OptionProcessor::format_options( std::ostream& aStream ) const {
	std::vector< std::string >::const_iterator
		group_name_i = m_option_group_names.begin(),
		group_name_end = m_option_group_names.end() ;
	for( ; group_name_i != group_name_end ; ++group_name_i ) {
		format_option_group( aStream, *group_name_i ) ;
	}
}

void OptionProcessor::format_option_group( std::ostream& aStream, std::string const& group_name ) const {
	aStream << group_name << ":\n" ;

	std::map< std::string, std::set< std::string > >::const_iterator
		group_i = m_option_groups.find( group_name ) ;
	assert( group_i != m_option_groups.end() ) ;
	std::set< std::string >::const_iterator
		option_i = group_i->second.begin(),
		end_option_i = group_i->second.end() ;

	std::size_t max_option_length = get_maximum_option_length() ;

	for( ; option_i != end_option_i; ++option_i ) {
		std::string const& option_name = *option_i ;
		aStream << " " ;
		format_option_and_description( aStream, option_name, max_option_length ) ;
	}

   	aStream << "\n";
}

std::string OptionProcessor::format_option_and_arguments( std::string const& option_name ) const {
	std::ostringstream oStream ;
	oStream
		<< option_name ;
	char arg_name = 'a' ;
	for( std::size_t i = 1 ; i <= (*this)[option_name].number_of_values_per_use(); ++i ) {
		oStream << " <" << arg_name++ << ">" ;
	}
	return oStream.str() ;
}

void OptionProcessor::format_option_and_description( std::ostream& aStream, std::string const& option_name, std::size_t max_option_length ) const {
	aStream
		<< std::setw( max_option_length+1 )
		<< std::right
		<< format_option_and_arguments( option_name ) ;

	aStream << ": " ;
	unsigned int current_column = max_option_length+4 ;

	aStream << wrap((*this)[option_name].description(), 100, current_column, current_column )
		<< "\n" ;
}

std::size_t OptionProcessor::get_maximum_option_length() const {
	OptionDefinitions::const_iterator
		option_i = m_option_definitions.begin(),
		end_option_i = m_option_definitions.end() ;
	std::size_t max_option_length = 0 ;
	for( ; option_i != end_option_i ; ++option_i ) {
		max_option_length = std::max( max_option_length, format_option_and_arguments( option_i->first ).size() ) ;
	}
	return max_option_length ;
}

std::size_t OptionProcessor::get_maximum_option_length( std::string const& group ) const {
	std::map< std::string, std::set< std::string > >::const_iterator
		group_i = m_option_groups.find( group ) ;
	assert( group_i != m_option_groups.end() ) ;
    std::set< std::string >::const_iterator
		option_i = group_i->second.begin(),
		end_option_i = group_i->second.end() ;
	std::size_t max_option_length = 0 ;
    for( ; option_i != end_option_i; ++option_i ) {
        max_option_length = std::max( max_option_length, format_option_and_arguments( *option_i ).size() ) ;
    }
	return max_option_length ;
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


