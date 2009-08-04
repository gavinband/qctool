#ifndef __GTOOL__OPTIONPROCESSOR_HPP__
#define __GTOOL__OPTIONPROCESSOR_HPP__


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

struct OptionProcessorHelpRequestedException: public std::exception {
	char const* what() const throw() { return "OptionProcessorHelpRequestedException" ; }
} ;

struct OptionProcessorMutuallyExclusiveOptionsSuppliedException: public std::exception {
	OptionProcessorMutuallyExclusiveOptionsSuppliedException( std::string const& option1, std::string const& option2 )
	: m_option1( option1 ),
	m_option2( option2 )
	{}
	
	~OptionProcessorMutuallyExclusiveOptionsSuppliedException() throw() {}

	char const* what() const throw() { return "OptionProcessorMutuallyExclusiveOptionsSuppliedException" ; }

	std::string const& first_option() const { return m_option1 ; }
	std::string const& second_option() const { return m_option2 ; }

private:
	
	std::string m_option1 ;
	std::string m_option2 ;
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

		void declare_group( std::string const& ) ;
		void set_help_option( std::string const& help_option_name ) { m_help_option_name = help_option_name ; }

		void option_excludes( std::string const& excluding_option, std::string const& excluded_option ) ;
		void option_excludes_group( std::string const& excluding_option, std::string const& excluded_option_group ) ;

		// Parse the options from argv, performing all needed checks.
		void process( int argc, char** argv ) ;

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

		friend std::ostream& operator<<( std::ostream& aStream, OptionProcessor::OptionDefinitions const& option_definitions ) ;
		friend std::ostream& operator<<( std::ostream& aStream, OptionProcessor const& options ) ;

	private:
		void calculate_option_groups() ;
		// Parse the options.  Store option values.  Ignore, but store unknown args for later reference
		void parse_options( int argc, char** argv ) ;
		bool try_to_parse_named_option_and_values( int argc, char** argv, int& i ) ;
		bool try_to_parse_positional_option( int argc, char** argv, int& i ) ;
		void process_unknown_options() ;
		void check_required_options_are_supplied() const ;
		void check_mutually_exclusive_options_are_not_supplied() const ;
		void preprocess_option_values() ;
		void check_option_values() ;
		std::string get_default_value( std::string const& arg ) const ;
		std::size_t get_maximum_option_length( std::string const& group ) const ;
		std::string format_option_and_arguments( std::string const& option_name ) const ;
		void format_options( std::ostream& ) const ;
		void format_option_group( std::ostream&, std::string const& ) const ;
		void format_option_and_description( std::ostream& aStream, std::string const& option_name, std::size_t max_option_length ) const ;
		
		OptionDefinitions m_option_definitions ;
		OptionValues m_option_values ;
		std::map< int, std::string > m_unknown_options ;
		std::string m_current_group ;
		std::map< std::string, std::set< std::string > > m_option_groups ;
		std::vector< std::string > m_option_group_names ;
		std::string m_help_option_name ;
		std::map< std::string, std::set< std::string > > m_option_exclusions ;
} ;


#endif

