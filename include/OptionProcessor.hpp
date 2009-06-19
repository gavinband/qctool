#ifndef __GTOOL__ARGUMENTPROCESSOR_HPP__
#define __GTOOL__ARGUMENTPROCESSOR_HPP__


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

struct ArgumentProcessingException: public GToolException
{
    ArgumentProcessingException( std::string msg )
    : GToolException( msg )
    {}
} ;

struct ArgumentParseException: public ArgumentProcessingException
{
    ArgumentParseException( std::string msg )
    : ArgumentProcessingException( msg )
    {}
} ;

struct ArgumentInvalidException: public ArgumentProcessingException
{
    ArgumentInvalidException( std::string msg )
    : ArgumentProcessingException( msg )
    {}
} ;

typedef std::map< std::string, ArgumentDefinition > ArgumentDefinitions ;

std::ostream& operator<<( std::ostream& aStream, ArgumentDefinitions const& argument_definitions ) {
    // Find longest argument and longest description
    ArgumentDefinitions::const_iterator
        defn_i = argument_definitions.begin(),
        defn_end = argument_definitions.end() ;

    std::size_t max_option_length = 0, max_description_length = 0 ;
     
    for( ; defn_i != defn_end; ++defn_i ) {
        max_option_length = std::max( max_option_length, defn_i->first.size()) ;
        max_description_length = std::max( max_description_length, defn_i->second.description().size() ) ;
    }

    defn_i = argument_definitions.begin() ;

    
    for( ; defn_i != defn_end; ++defn_i ) {
        aStream << std::setw( max_option_length+1 )
                << std::right
                << defn_i->first
                << " " ;
				
		std::string flags = "" ;
        if( defn_i->second.is_required() ) {
            flags += "required" ;
		}
        if( defn_i->second.takes_value() ) {
			if( flags.size() > 0 ) {
				flags += "," ;
			}
			flags += "takes value" ;
		}
        aStream << ":\t"
                << std::setw( max_description_length+1 )
                << std::left
                << defn_i->second.description()
                << "\t"
				<< "(" + flags + ")"
				<< ".\n" ;
    }

    return aStream ;
}


struct OptionProcessor {
	public:
		OptionProcessor() {}
		OptionProcessor( OptionProcessor const& other )
			: m_argument_definitions( other.m_argument_definitions ),
				m_argument_values( other.m_argument_values ),
				m_unknown_arguments( other.m_unknown_arguments )
		{}

		~OptionProcessor() {} ;
		
		ArgumentDefinitions const& argument_definitions() const { return m_argument_definitions ; }

		ArgumentDefinition& operator[]( std::string const& arg ) {
			return m_argument_definitions[ arg ] ;
		}

		ArgumentDefinition const& operator[] ( std::string const& arg ) const {
			ArgumentDefinitions::const_iterator defn_i = m_argument_definitions.find( arg ) ;
			if( defn_i == m_argument_definitions.end() )
				throw ArgumentInvalidException( "Option \"" + arg + "\" is not defined." ) ;
			return defn_i->second ;
		}

		// Parse the options from argv, performing all needed checks.
		void process( int argc, char** argv ) {
			parse_arguments( argc, argv ) ;
			check_required_arguments_are_supplied() ;
			check_argument_values() ;
		}

		// Parse the arguments.  Store argument values.  Ignore, but store unknown args for later reference
		void parse_arguments( int argc, char** argv ) {
			int i = 0;
			while( ++i < argc ) {
				std::string arg = argv[i] ;
				std::map< std::string, ArgumentDefinition >::const_iterator arg_i = m_argument_definitions.find( arg ) ;
				if( arg_i != m_argument_definitions.end() ) {
					if( arg_i->second.takes_value() ) {
						if( ++i == argc ) {
							throw ArgumentParseException( "Option \"" + arg + "\" requires a value." ) ;
						}
						else {
							m_argument_values[ arg ] = argv[i] ;
						}
					}
					else {
						m_argument_values[ arg ] = "1" ;
					}
				}
				else
				{
					m_unknown_arguments.insert( argv[i] ) ;
				}
			}
			
			if( !m_unknown_arguments.empty()) {
				std::string unknown_args ;
				std::set< std::string >::const_iterator i = m_unknown_arguments.begin() ;
				for( ; i != m_unknown_arguments.end(); ++i ) {
					if( i != m_unknown_arguments.begin() )
						unknown_args += ", " ;
					unknown_args += *i ;
				}
				throw ArgumentParseException( "The following options were not recognised: " + unknown_args + "." ) ;
			}
		}

		void check_required_arguments_are_supplied() {
			std::map< std::string, ArgumentDefinition >::const_iterator
				defn_i = m_argument_definitions.begin(), 
				defn_end = m_argument_definitions.end() ;

			for( ; defn_i != defn_end; ++defn_i ) {
				if( defn_i->second.is_required() && !check_if_argument_was_supplied( defn_i->first ) && !defn_i->second.has_default_value()) {
					throw ArgumentProcessingException( "Option \"" + defn_i->first + "\" must be supplied." ) ;
				}
			}
		}

		void check_argument_values() {
			std::map< std::string, std::string >::const_iterator
				arg_i = m_argument_values.begin(),
				arg_end = m_argument_values.end() ;

			for( ; arg_i != arg_end; ++arg_i ) {
				ArgumentDefinitions::const_iterator defn_i
					= m_argument_definitions.find( arg_i->first ) ;
				if( defn_i->second.value_checker() != 0 ) {
				   defn_i->second.value_checker()( arg_i->first, arg_i->second ) ; 
				}
			}
		}

		// check if the given argument (which must be valid) was supplied.
		bool check_if_argument_was_supplied( std::string const& arg ) const {
			return ( m_argument_values.find( arg ) != m_argument_values.end() ) ;
		}

		// get the value of the given argument as the given type.
		template< typename T >
		T get_value( std::string const& arg ) const {
			std::istringstream s( get_value< std::string >( arg )) ;
			T t ;
			s >> t ;
			return t ;
		}

		std::string get_default_value( std::string const& arg ) const {
			std::map< std::string, ArgumentDefinition >::const_iterator defn_i
				= m_argument_definitions.find( arg ) ;
			assert( defn_i != m_argument_definitions.end() ) ;
			assert( defn_i->second.has_default_value() ) ;
			return defn_i->second.default_value() ;
		}
	
    private:

		std::map< std::string, ArgumentDefinition > m_argument_definitions ;
		std::map< std::string, std::string > m_argument_values ;
		std::set< std::string > m_unknown_arguments ;
} ;

std::ostream& operator<<( std::ostream& aStream, OptionProcessor const& options ) {
	return aStream << "Options:\n" << options.argument_definitions() ;
}

// get the value of the given argument as a string.
// We specialise this for strings to handle the case where
// the value contains spaces.
template<>
std::string OptionProcessor::get_value( std::string const& arg ) const {
	std::map< std::string, std::string >::const_iterator arg_i ;
	arg_i = m_argument_values.find( arg ) ;
	
	if( arg_i == m_argument_values.end() ) {
		return get_default_value( arg ) ;
	}

	return arg_i->second ;
}

#endif

