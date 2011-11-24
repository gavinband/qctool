#ifndef GENFILE_SNP_IDENTIFYING_DATA_HPP
#define GENFILE_SNP_IDENTIFYING_DATA_HPP

#include <string>
#include <vector>
#include "genfile/GenomePosition.hpp"

namespace genfile {
	struct SNPIdentifyingData
	{
	public:
		SNPIdentifyingData() ;
		SNPIdentifyingData(
			std::string const& SNPID,
			std::string const& RSID,
			GenomePosition const& position,
			char first_allele,
			char second_allele
		) ;
		SNPIdentifyingData(
			std::string const& SNPID,
			std::string const& RSID,
			GenomePosition const& position,
			std::string const& first_allele,
			std::string const& second_allele
		) ;
		SNPIdentifyingData( SNPIdentifyingData const& other ) ;
		SNPIdentifyingData& operator=( SNPIdentifyingData const& other ) ;

		std::string& SNPID() { return m_SNPID ; }
		std::string& rsid() { return m_RSID ;}
		GenomePosition& position() { return m_position ;}
		std::string& first_allele() { return m_first_allele ;}
		std::string& second_allele() { return m_second_allele ;}

		std::string const& get_SNPID() const { return m_SNPID ;}
		std::string const& get_rsid() const { return m_RSID ;}
		GenomePosition const& get_position() const { return m_position ;}
		std::string const& get_first_allele() const { return m_first_allele ;}
		std::string const& get_second_allele() const { return m_second_allele ;}
	public:
		struct CompareFields {
			CompareFields( std::string const& fields_to_compare ) ;
			CompareFields( CompareFields const& other ) ;

			enum { eSNPID = 0x1, eRSID = 0x2, ePosition = 0x4, eAlleles = 0x8, eMask = 0xF } ;
			bool operator()( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) const ;
			bool are_equal( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) const ;
			bool check_if_comparable_fields_are_known( SNPIdentifyingData const& value ) const ;

		private:
			static std::vector< int > parse_fields_to_compare( std::string const& field_spec ) ;
			
			std::vector< int > const m_fields_to_compare ;
		} ;
	private:
		std::string m_SNPID ;
		std::string m_RSID ;
		GenomePosition m_position ;
		std::string m_first_allele ;
		std::string m_second_allele ;
		
		CompareFields& operator=( CompareFields const& other ) ;
	} ;	
	
	std::ostream& operator<<( std::ostream&, SNPIdentifyingData const& ) ;
	std::ostream& operator<<( std::ostream&, std::vector< SNPIdentifyingData > const& ) ;
	bool operator==( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) ;
	bool operator!=( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) ;
	bool operator<( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) ;
}

#endif
