#ifndef __GTOOL__OPTION_DEFINITION_HPP__
#define __GTOOL__OPTION_DEFINITION_HPP__


#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <cassert>

struct OptionDefinition {

	public:
		OptionDefinition() ;
		OptionDefinition( OptionDefinition const& other ) ;

		typedef void (*value_checker_t)( std::string const&, std::vector< std::string > const& ) ;
		typedef std::vector< std::string > (*value_preprocessor_t)( std::string const&, std::vector< std::string > const& ) ;

		std::string const& group() const { return m_group ; }
		std::string description() const { return m_description ; } 
		bool is_required() const { return m_is_required ; }
		unsigned int number_of_values_per_use() const { return m_number_of_values_per_use ; }
		unsigned int maximum_number_of_repeats() const {
			return ( m_number_of_values_per_use == 0 ) ? 0 : ( m_upper_number_of_permitted_values / m_number_of_values_per_use ) ;
		}
		
		unsigned int minimum_number_of_permitted_values() const { return m_upper_number_of_permitted_values ; }
		unsigned int maximum_number_of_permitted_values() const { return m_upper_number_of_permitted_values ; }
		bool takes_values() const { return m_upper_number_of_permitted_values > 0 ; }
		std::vector< value_checker_t > value_checkers() const { return m_value_checkers ; } 
		std::vector< value_preprocessor_t > value_preprocessors() const { return m_value_preprocessors ; } 
		bool has_default_value() const { return m_has_default_value ; }
		std::string default_value() const { return m_default_value ; }
		bool takes_value_by_position() const { return m_position > 0 ; }
		int position() const { return m_position ; }
		bool hidden() const { return m_hidden ; }
		
		OptionDefinition& set_group( std::string const& group ) { m_group = group ; return *this ; }
		OptionDefinition& set_description( std::string desc ) { m_description = desc ; return *this ; }
		OptionDefinition& set_is_required() { m_is_required = true ; return *this ; }
		OptionDefinition& set_takes_single_value() {
			m_number_of_values_per_use = 1u ;
			m_lower_number_of_permitted_values = 1u ;
			m_upper_number_of_permitted_values = 1u ;
			return *this ;
		}
		OptionDefinition& set_takes_values() {
			m_number_of_values_per_use = std::max( m_number_of_values_per_use, 1u ) ;
			m_lower_number_of_permitted_values = std::max( m_lower_number_of_permitted_values, 1u ) ;
			m_upper_number_of_permitted_values = std::max( m_upper_number_of_permitted_values, 1u ) ;
			return *this ;
		}
		OptionDefinition& set_number_of_values_per_use( unsigned int n ) {
			m_number_of_values_per_use = n;
			m_upper_number_of_permitted_values = std::max( m_upper_number_of_permitted_values, n ) ;
			return *this ;
		}
		OptionDefinition& set_maximum_number_of_repeats( unsigned int n ) {
			m_upper_number_of_permitted_values = n * m_number_of_values_per_use ;
			return *this ;
		}
		OptionDefinition& add_value_checker( value_checker_t value_checker ) {
			assert( value_checker != 0 ) ;
			m_value_checkers.push_back( value_checker ) ; 
			return *this ; 
		}
		OptionDefinition& add_value_preprocessor( value_preprocessor_t value_preprocessor ) {
			assert( value_preprocessor != 0 ) ;
			m_value_preprocessors.push_back( value_preprocessor );
			return *this ;
		}
		template< typename T > OptionDefinition& set_default_value( T const& value ) {
			std::ostringstream aStream ;
			aStream << value ;
			m_default_value = aStream.str() ;
			m_has_default_value = true ;
			return *this ;
		}
		OptionDefinition& set_takes_value_by_position( int position ) {
			m_position = position ;
			return *this ;
		}
		OptionDefinition& set_hidden() {
			m_hidden = true ;
		}

		std::vector< std::string > preprocess_option_values( std::string const& option_name, std::vector< std::string > const& option_values ) const ;
		void check_option_values( std::string const& option_name, std::vector< std::string > const& option_values ) const ;

    private:

		std::string m_group ;
		std::string m_description ;
		bool m_is_required ;
		unsigned int m_number_of_values_per_use, m_lower_number_of_permitted_values, m_upper_number_of_permitted_values ;
		std::vector< value_checker_t > m_value_checkers ;
		std::vector< value_preprocessor_t > m_value_preprocessors ;
		bool m_has_default_value ;
		std::string m_default_value ;
		int m_position ;
		bool m_hidden ;
} ;

typedef OptionDefinition OptionDefinition ;

#endif

