#ifndef __GTOOL__OPTION_DEFINITION_HPP__
#define __GTOOL__OPTION_DEFINITION_HPP__


#include <string>
#include <sstream>
#include <cassert>

struct ArgumentDefinition {

	public:
		ArgumentDefinition()
			: m_description(""),
				m_is_required( false ),
				m_takes_value( false ),
				m_value_checker(0),
				m_has_default_value( false ),
				m_default_value( "" )
		{}

		ArgumentDefinition( ArgumentDefinition const& other )
			: m_description( other.m_description ),
				m_is_required( other.m_is_required ),
				m_takes_value( other.m_takes_value ),
				m_value_checker( other.m_value_checker ),
				m_has_default_value( other.m_has_default_value ),
				m_default_value( other.m_default_value )
		{}
		
		typedef void (*value_checker_t)( std::string const&, std::string const& ) ;
	 
		ArgumentDefinition& set_description( char const* desc ) { m_description = desc ; return *this ; }
		ArgumentDefinition& set_is_required() { m_is_required = true ; return *this ; }
		ArgumentDefinition& set_takes_value() { m_takes_value = true ; return *this ; }
		ArgumentDefinition& set_value_checker( value_checker_t value_checker ) { m_value_checker = value_checker ; return *this ; }
		template< typename T > ArgumentDefinition& set_default_value( T const& value ) {
			std::ostringstream aStream ;
			aStream << value ;
			m_default_value = aStream.str() ;
			m_has_default_value = true ;
			m_takes_value = true ;
			return *this ;
		}
		
		std::string description() const { return m_description ; } 
		bool is_required() const { return m_is_required ; }
		bool takes_value() const { return m_takes_value ; }
		value_checker_t value_checker() const { return m_value_checker ; } 
		bool has_default_value() const { return m_has_default_value ; }
		std::string default_value() const { return m_default_value ; }

    private:

		char const* m_description ;
		bool m_is_required ;
		bool m_takes_value ;
		value_checker_t m_value_checker ;
		bool m_has_default_value ;
		std::string m_default_value ;
} ;


#endif

