#include <string>
#include <sstream>
#include <cassert>
#include "appcontext/OptionDefinition.hpp"
#include "appcontext/OptionProcessor.hpp"

namespace appcontext {
	OptionDefinition::OptionDefinition( std::string const& name ):
		m_name( name ),
		m_number_of_values_per_use( 0 ),
		m_minimum_multiplicity( 0 ),
		m_maximum_multiplicity( 1 ),
		m_default_values(),
		m_position( -1 ),
		m_hidden( false )
	{}

	OptionDefinition::OptionDefinition( OptionDefinition const& other ):
		m_name( other.m_name ),
		m_group( other.m_group ),
		m_description( other.m_description ),
		m_number_of_values_per_use( other.m_number_of_values_per_use ),
		m_minimum_multiplicity( other.m_minimum_multiplicity ),
		m_maximum_multiplicity( other.m_maximum_multiplicity ),
		m_default_values( other.m_default_values ),
		m_position( other.m_position ),
		m_hidden( false )
	{}

	OptionDefinition& OptionDefinition::operator=( OptionDefinition const& other ) {
		m_name = other.m_name ;
		m_group = other.m_group ;
		m_description = other.m_description ;
		m_number_of_values_per_use = other.m_number_of_values_per_use ;
		m_minimum_multiplicity = other.m_minimum_multiplicity ;
		m_maximum_multiplicity = other.m_maximum_multiplicity ;
		m_default_values = other.m_default_values ;
		m_position = other.m_position ;
		m_hidden = other.m_hidden ;
		return *this ;
	}

	void OptionDefinition::check_option_values( std::string const& option_name, std::vector< std::string > const& option_values ) const {
		// Check if number or values is too few...
		if( number_of_values_per_use() != eUntilNextOption ) {
			if( option_values.size() < m_minimum_multiplicity * number_of_values_per_use() ) {
				std::ostringstream ostr ;
				ostr << "Option \"" << option_name << "\" requires at least " << m_minimum_multiplicity * number_of_values_per_use() << " value" ;
				if( m_minimum_multiplicity * number_of_values_per_use() > 1 ) ostr << "s" ;
		
				throw OptionValueInvalidException( option_name, option_values, ostr.str() ) ;
			}

			// ...or too many.
			if( option_values.size() > m_maximum_multiplicity * number_of_values_per_use() ) {
				std::ostringstream ostr ;
				ostr << "Option \"" << option_name << "\" takes at most " << m_maximum_multiplicity * number_of_values_per_use() << " value" ;
				if( m_maximum_multiplicity * number_of_values_per_use() > 1 ) ostr << "s" ;

				throw OptionValueInvalidException( option_name, option_values, ostr.str() ) ;
			}
		}

		std::vector< value_checker_t >::const_iterator
			i = m_value_checkers.begin(),
			end_i = m_value_checkers.end() ;
		for( ; i != end_i; ++i ) {
			(*i)( option_name, option_values ) ;
		}
	}
	
	std::string OptionDefinition::format_spec( std::size_t indent ) const {
		std::ostringstream ostream ;
		std::string const space( indent, ' ' ) ;
		std::string const space2 = space + space ;
		ostream << space << "\"" << m_name << "\": {\n" ;
		ostream
			<< space2 << "\"group\": \"" << m_group << "\",\n"
			<< space2 << "\"description\": \"" << m_description << "\",\n"
			<< space2 << "\"is_required\": " << m_is_required << ",\n"
			<< space2 << "\"multiplicity\": {\n"
			<< space2 << space << "\"min\": " << m_minimum_multiplicity << ",\n"
			<< space2 << space << "\"max\": " << m_maximum_multiplicity << "\n"
			<< space2 << "},\n" ;
		if( m_number_of_values_per_use == eUntilNextOption ) {
			ostream << space2 << "\"values_per_use\": \"untilNextOption\",\n" ;
		} else {
			ostream << space2 << "\"values_per_use\": " << m_number_of_values_per_use << ",\n" ;
		}
		ostream << space2 << "\"default\": [" ;
		for( std::size_t i = 0; i < m_default_values.size(); ++i ) {
			ostream << ((i>0) ? ", " : " " ) << "\"" << m_default_values[i] << "\"" ;
		}
		ostream
			<< ((m_default_values.size() > 0) ? " ]\n" : "]\n" )
			<< space << "}" ;
		return ostream.str() ;
	}
}
