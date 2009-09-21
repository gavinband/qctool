#include <set>
#include <cassert>
#include <fstream>
#include <iostream>

#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SNPIDMatchesCondition.hpp"
#include "GenRowStatistics.hpp"
#include "parse_utils.hpp"

SNPIDMatchesCondition::SNPIDMatchesCondition( std::string const& expression )
 : m_expression( expression ),
	m_wildcard_char( '*' )
{
	setup() ;
}

void SNPIDMatchesCondition::setup() {
	std::size_t pos = m_expression.find( m_wildcard_char ) ;
	if( pos == std::string::npos ) {
		m_prefix = m_expression ;
		m_suffix = "" ;
	}
	else {
		m_prefix = m_expression.substr( 0, pos ) ;
		m_suffix = m_expression.substr( pos + 1, m_expression.size() ) ;
	}
}

bool SNPIDMatchesCondition::check_if_satisfied( string_to_value_map const& statistics ) const {
	GenRowStatistics const* row_statistics_ptr = dynamic_cast< GenRowStatistics const* >( &statistics ) ;
	if( !row_statistics_ptr ) {
		throw ConditionException( "SNPIDMatchesCondition only supports GenRowStatistics." ) ;
	}
	std::string const& SNPID = row_statistics_ptr->row().SNPID() ;
	return string_has_prefix_and_suffix( SNPID, m_prefix, m_suffix ) ;
}

void SNPIDMatchesCondition::format_to_stream( std::ostream& oStream ) const {
	oStream << "SNPID-matches(" << m_expression << ")" ;
}