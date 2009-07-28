#ifndef __GTOOL__OPTIONPROCESSOR_HPP__
#define __GTOOL__OPTIONPROCESSOR_HPP__


#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "GenRow.hpp"
#include "GToolException.hpp"
#include "OptionDefinition.hpp"

struct OptionProcessingException: public GToolException
{
	OptionProcessingException( std::string option, std::vector< std::string > values, std::string msg ) ;
	OptionProcessingException( std::string option, std::string msg ) ;
	~OptionProcessingException() throw() ;
	
	char const* what() const throw() ;
private:
	
	std::string m_option ;
	std::vector< std::string > m_values ;
	std::string m_msg ;
} ;

struct OptionParseException: public OptionProcessingException
{
	OptionParseException( std::string option, std::vector< std::string > values, std::string msg ) ;
} ;

struct OptionValueInvalidException: public OptionProcessingException
{
	OptionValueInvalidException( std::string option, std::vector< std::string > values, std::string msg ) ;
} ;

class OptionProcessor {
		typedef std::map< std::string, OptionDefinition > OptionDefinitions ;
		typedef std::map< std::string, std::vector< std::string > > OptionValues ; 

	public:
		OptionProcessor() ;
		OptionProcessor( OptionProcessor const& other ) ;
		~OptionProcessor() ;
		
		OptionDefinitions const& option_definitions() const ;
		OptionDefinition& operator[]( std::string const& arg ) ;
		OptionDefinition const& operator[] ( std::string const& arg ) const ;

		// Parse the options from argv, performing all needed checks.
		void process( int argc, char** argv ) ;

		// Parse the options.  Store option values.  Ignore, but store unknown args for later reference
		void parse_options( int argc, char** argv ) ;
		bool try_to_parse_named_option_and_values( int argc, char** argv, int& i ) ;
		bool try_to_parse_positional_option( int argc, char** argv, int& i ) ;
		void process_unknown_options() ;
		void check_required_options_are_supplied() ;
		void preprocess_option_values() ;
		void check_option_values() ;

		// check if the given option (which must be valid) was supplied.
		bool check_if_option_was_supplied( std::string const& arg ) const ;

		// get the value of the given option as the given type.
		template< typename T >
		T get_value( std::string const& arg ) const {
			std::istringstream s( get_value< std::string >( arg )) ;
			T t ;
			s >> t ;
			return t ;
		}

		// get the value of the given option as the given type.
		template< typename T >
		std::vector< T > get_values( std::string const& arg ) const {
			std::vector< std::string > values = get_values< std::string >( arg ) ;
			std::vector<T> result ;
			result.reserve( values.size() ) ;
			for( std::size_t i = 0; i < values.size(); ++i ) {
				std::istringstream istr( values[i] ) ;
				T t ;
				istr >> t ;
				result.push_back( t ) ;
			}
			return result ;
		}

		std::string get_default_value( std::string const& arg ) const ;
	
    private:
		OptionDefinitions m_option_definitions ;
		OptionValues m_option_values ;
		std::map< int, std::string > m_unknown_options ;

	public:
		friend std::ostream& operator<<( std::ostream& aStream, OptionProcessor::OptionDefinitions const& option_definitions ) ;
		friend std::ostream& operator<<( std::ostream& aStream, OptionProcessor const& options ) ;
} ;


#endif

