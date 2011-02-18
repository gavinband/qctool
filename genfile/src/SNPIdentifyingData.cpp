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

	std::ostream& operator<<( std::ostream& out, std::vector< SNPIdentifyingData > const& data ) {
		for( std::size_t i = 0; i < data.size(); ++i ) {
			out << data[i] << "\n" ;
		}
		return out ;
	}
	
	bool operator==( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) {
		return left.get_SNPID() == right.get_SNPID() &&
			left.get_rsid() == right.get_rsid() &&
			left.get_position() == right.get_position() &&
			left.get_first_allele() == right.get_first_allele() &&
			left.get_second_allele() == right.get_second_allele() ;
	}

	bool operator!=( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) {
		return left.get_SNPID() != right.get_SNPID() ||
			left.get_rsid() != right.get_rsid() ||
			left.get_position() != right.get_position() ||
			left.get_first_allele() != right.get_first_allele() ||
			left.get_second_allele() != right.get_second_allele() ;
	}
	
    bool operator<( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) {
		return(
			(left.get_position() < right.get_position())
			||
			(
				(left.get_position() == right.get_position())
				&&
				(
					( left.get_rsid() < right.get_rsid() )
					||
					(
						(left.get_rsid() == right.get_rsid())
						&&
						(
							(left.get_SNPID() < right.get_SNPID())
							||
							(
								(left.get_SNPID() == right.get_SNPID())
								&&
								(
									(left.get_first_allele() < right.get_first_allele())
									||
									(
										(left.get_first_allele() == right.get_first_allele())
										&&
										(left.get_second_allele() < right.get_second_allele())
									)
								)
							)
						)
					)
				)
			)
		) ;
	}
}
