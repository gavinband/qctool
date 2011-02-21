#include <string>
#include <vector>
#include <cassert>
#include "genfile/GenomePosition.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/string_utils.hpp"

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
	
	SNPIdentifyingData::CompareFields::CompareFields( std::string const& fields_to_compare ):
		m_fields_to_compare( parse_fields_to_compare( fields_to_compare ) )
	{
		assert( !m_fields_to_compare.empty() ) ;
	}
	
	SNPIdentifyingData::CompareFields::CompareFields( CompareFields const& other ):
		m_fields_to_compare( other.m_fields_to_compare )
	{}

	std::vector< int > SNPIdentifyingData::CompareFields::parse_fields_to_compare( std::string const& field_spec ) {
		std::vector< std::string > elts = string_utils::split_and_strip( field_spec, ",", " \t\n\r" ) ;
		assert( elts.size() > 0 ) ;
		std::vector< int > result ;
		for( std::size_t i = 0; i < elts.size(); ++i ) {
			if( elts[i] == "SNPID" ) {
				result.push_back( eSNPID ) ;
			}
			else if( elts[i] == "rsid" ) {
				result.push_back( eRSID ) ;
			}
			else if( elts[i] == "position" ) {
				result.push_back( ePosition ) ;
			}
			else if( elts[i] == "alleles" ) {
				result.push_back( eAlleles ) ;
			}
			else {
				assert(0) ;
			}
		}
		return result ;
	}

	bool SNPIdentifyingData::CompareFields::operator()( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) const {
		// Lexicographical compare in the same order as operator< above, but using only the specified fields.

		for( std::size_t i = 0; i < m_fields_to_compare.size(); ++i ) {
			switch( m_fields_to_compare[i] ) {
				case ePosition:
					if( left.get_position() > right.get_position() ) {
						return false ;
					}
					else if( left.get_position() < right.get_position() ) {
						return true ;
					}
					break ;
				case eRSID:
					if( left.get_rsid() > right.get_rsid() ) {
						return false ;
					}
					else if( left.get_rsid() < right.get_rsid() ) {
						return true ;
					}
					break ;
				case eSNPID:
					if( left.get_SNPID() > right.get_SNPID() ) {
						return false ;
					}
					else if( left.get_SNPID() < right.get_SNPID() ) {
						return true ;
					}
					break ;
				case eAlleles:
					if( left.get_first_allele() > right.get_first_allele() ) {
						return false ;
					}
					else if( left.get_first_allele() < right.get_first_allele() ) {
						return true ;
					}
					else if( left.get_second_allele() > right.get_second_allele() ) {
						return false ;
					}
					else if( left.get_second_allele() < right.get_second_allele() ) {
						return true ;
					}
					break ;
				default:
					assert(0) ;
					break ;
			}
		}
		return false ;
	}
	
	bool SNPIdentifyingData::CompareFields::are_equal( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) const {
		for( std::size_t i = 0; i < m_fields_to_compare.size(); ++i ) {
			switch( m_fields_to_compare[i] ) {
				case ePosition:
					if( left.get_position() != right.get_position() ) {
						return false ;
					}
					break ;
				case eRSID:
					if( left.get_rsid() != right.get_rsid() ) {
						return false ;
					}
					break ;
				case eSNPID:
					if( left.get_SNPID() != right.get_SNPID() ) {
						return false ;
					}
					break ;
				case eAlleles:
					if( left.get_first_allele() != right.get_first_allele() ) {
						return false ;
					}
					if( left.get_second_allele() != right.get_second_allele() ) {
						return false ;
					}
				default:
					assert(0) ;
					break ;
			}
		}
		return true ;
	}
	
}
