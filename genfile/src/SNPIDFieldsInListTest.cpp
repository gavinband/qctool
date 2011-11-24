#include <set>
#include <string>
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPIDFieldsInListTest.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	SNPIDFieldsInListTest::SNPIDFieldsInListTest( std::set< std::string > id_fields )
	: m_id_fields( id_fields )
	{}

	bool SNPIDFieldsInListTest::operator()(
		std::string SNPID,
		std::string RSID,
		GenomePosition,
		std::string,
		std::string
	) const {
		return (m_id_fields.find( SNPID ) != m_id_fields.end()) || (m_id_fields.find( RSID ) != m_id_fields.end()) ;
	}
	
	std::string SNPIDFieldsInListTest::display() const {
		std::string result = "SNPID or RSID in { " ;
		if( m_id_fields.size() <= 10 ) {
			for( std::set< std::string >::const_iterator i = m_id_fields.begin(); i != m_id_fields.end(); ++i ) {
				if( i != m_id_fields.begin() ) {
					result += ", " ;
				}
				result += *i ;
			}
		}
		else {
			result += "set of " + string_utils::to_string( m_id_fields.size() ) ;
		}
		return result + " }";
	}
}
