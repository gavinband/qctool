#include <set>
#include <string>
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPIDFieldsInListTest.hpp"

namespace genfile {
	SNPIDFieldsInListTest::SNPIDFieldsInListTest( std::set< std::string > id_fields )
	: m_id_fields( id_fields )
	{}

	bool SNPIDFieldsInListTest::operator()(
		std::string SNPID,
		std::string RSID,
		GenomePosition,
		char,
		char
	) const {
		std::string const key = SNPID + ":" + RSID ;
		return m_id_fields.find( key ) != m_id_fields.end() ;
	}
}
