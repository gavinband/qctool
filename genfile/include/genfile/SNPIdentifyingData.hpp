#ifndef GENFILE_SNP_IDENTIFYING_DATA_HPP
#define GENFILE_SNP_IDENTIFYING_DATA_HPP

#include <string>
#include <vector>
#include "genfile/GenomePosition.hpp"

namespace genfile {
	struct SNPIdentifyingData
	{
		SNPIdentifyingData() ;
		SNPIdentifyingData( std::string const& SNPID, std::string const& RSID, GenomePosition const& position, char first_allele, char second_allele ) ;
		SNPIdentifyingData( SNPIdentifyingData const& other ) ;
		SNPIdentifyingData& operator=( SNPIdentifyingData const& other ) ;

		std::string& SNPID() { return m_SNPID ; }
		std::string& rsid() { return m_RSID ;}
		GenomePosition& position() { return m_position ;}
		char& first_allele() { return m_first_allele ;}
		char& second_allele() { return m_second_allele ;}

		std::string const& get_SNPID() const { return m_SNPID ;}
		std::string const& get_rsid() const { return m_RSID ;}
		GenomePosition const& get_position() const { return m_position ;}
		char get_first_allele() const { return m_first_allele ;}
		char get_second_allele() const { return m_second_allele ;}

	private:
		std::string m_SNPID ;
		std::string m_RSID ;
		GenomePosition m_position ;
		char m_first_allele ;
		char m_second_allele ;
	} ;	
	
	std::ostream& operator<<( std::ostream&, SNPIdentifyingData const& ) ;
	std::ostream& operator<<( std::ostream&, std::vector< SNPIdentifyingData > const& ) ;
	bool operator==( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) ;
	bool operator!=( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) ;
	bool operator<( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) ;
}

#endif
