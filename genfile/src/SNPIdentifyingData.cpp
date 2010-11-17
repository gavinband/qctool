#include <string>
#include "genfile/GenomePosition.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	SNPIdentifyingData::SNPIdentifyingData() {}
	
	SNPIdentifyingData::SNPIdentifyingData(
		std::string const& SNPID,
		std::string const& RSID,
		GenomePosition const& position,
		char first_allele,
		char second_allele
	):
		m_SNPID( SNPID ),
		m_RSID( RSID ),
		m_position( position ),
		m_first_allele( first_allele ),
		m_second_allele( second_allele )
	{}
	
	SNPIdentifyingData::SNPIdentifyingData( SNPIdentifyingData const& other ):
		m_SNPID( other.m_SNPID ),
		m_RSID( other.m_RSID ),
		m_position( other.m_position ),
		m_first_allele( other.m_first_allele ),
		m_second_allele( other.m_second_allele )
	{}

	SNPIdentifyingData& SNPIdentifyingData::operator=( SNPIdentifyingData const& other ) {
		m_SNPID = other.m_SNPID ;
		m_RSID = other.m_RSID ;
		m_position = other.m_position ;
		m_first_allele = other.m_first_allele ;
		m_second_allele = other.m_second_allele ;
		return *this ;
	}
	
	std::ostream& operator<<( std::ostream& out, SNPIdentifyingData const& data ) {
		return out << data.get_SNPID()
			<< " " << data.get_rsid()
			<< " " << data.get_position()
			<< " " << data.get_first_allele()
			<< " " << data.get_second_allele() ;
	}
	
}
