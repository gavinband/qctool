#include <string>
#include <cassert>
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPIDMatchesTest.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	SNPIDMatchesTest::SNPIDMatchesTest( std::string const& expression ) {
		setup( expression ) ;
	}

	void SNPIDMatchesTest::setup( std::string const& expression ) {
		std::string lowered_expression = string_utils::to_lower( expression ) ;
		std::vector< std::string > bits = string_utils::split_and_strip( lowered_expression, "~" ) ;
		if( bits.size() == 1 ) {
			bits.insert( bits.begin(), "either" ) ;
		}
		if( bits.size() != 2 || ( bits[0] != "snpid" && bits[0] != "rsid" && bits[0] != "either" )) {
			throw BadArgumentError( "SNPIDMatchesTest::setup()", "expression = \"" + expression + "\"" ) ;
		}
		if( bits[0] == "snpid" ) {
			m_type = eSNPID ;
		}
		else if( bits[0] == "rsid" ) {
			m_type = eRSID ;
		}
		else if( bits[0] == "either" ) {
			m_type = eEITHER ;
		}

		m_expression = bits[1] ;
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

	bool SNPIDMatchesTest::operator()(
		std::string SNPID,
		std::string RSID,
		GenomePosition,
		char,
		char
	) const {
		if( m_type == eSNPID ) {
			return match( SNPID ) ;
		}
		else if( m_type == eRSID ) {
			return match( RSID ) ;
		}
		else if( m_type == eEITHER ) {
			return match( SNPID ) || match( RSID ) ;
		}
		else {
			assert(0) ;
		}
	}
	
	bool SNPIDMatchesTest::match( std::string const& s ) const {
		return s.size() >= (m_prefix.size() + m_suffix.size()) && s.compare( 0, m_prefix.size(), m_prefix ) == 0 && s.compare( s.size() - m_suffix.size(), m_suffix.size(), m_suffix ) == 0 ;
	}
	
	std::string SNPIDMatchesTest::display() const {
		if( m_type == eSNPID ) {
			return "SNPID~" + m_expression ;
		}
		else if( m_type == eRSID ) {
			return "RSID~" + m_expression ;
		}
		else if( m_type == eEITHER ) {
			return "(SNPID or RSID)~" + m_expression ;
		}
		else {
			assert(0) ;
		}
	}
}
