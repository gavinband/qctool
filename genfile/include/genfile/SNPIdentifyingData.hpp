
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNP_IDENTIFYING_DATA_HPP
#define GENFILE_SNP_IDENTIFYING_DATA_HPP

#include <string>
#include <vector>
#include "genfile/GenomePosition.hpp"
#include "genfile/string_utils/slice.hpp"

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
		void swap_alleles() ;
		std::size_t number_of_alleles() const { return 2 ; }

		void set_SNPID( std::string const& SNPID ) { m_SNPID = SNPID ; }
		void set_rsid( std::string const& rsid ) { m_RSID = rsid ;}
		void set_position( GenomePosition const& position ) { m_position = position ;}
		void set_first_allele( std::string const& allele ) { m_first_allele = allele ;}
		void set_second_allele( std::string const& allele ) { m_second_allele = allele ;}

		std::string const& get_SNPID() const { return m_SNPID ;}
		std::string const& get_rsid() const { return m_RSID ;}
		GenomePosition const& get_position() const { return m_position ;}
		std::string const& get_first_allele() const { return m_first_allele ;}
		std::string const& get_second_allele() const { return m_second_allele ;}
	public:
		struct CompareFields {
			CompareFields() ;
			CompareFields( std::string const& fields_to_compare, bool flip_alleleles_if_necessary = false ) ;
			CompareFields( CompareFields const& other ) ;

			CompareFields& operator=( CompareFields const& other ) ;

			bool get_flip_alleles_if_necessary() const { return m_flip_alleles_if_necessary ; }
			std::vector< int > const& get_compared_fields() const { return m_fields_to_compare ;}
			
			enum { eSNPID = 0x1, eRSID = 0x2, ePosition = 0x4, eAlleles = 0x8, eMask = 0xF } ;
			bool operator()( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) const ;
			bool are_equal( SNPIdentifyingData const& left, SNPIdentifyingData const& right ) const ;
			bool check_if_comparable_fields_are_known( SNPIdentifyingData const& value ) const ;

			std::string get_summary() const ;

		private:
			static std::vector< int > parse_fields_to_compare( std::string const& field_spec ) ;
			std::vector< int > m_fields_to_compare ;
			bool m_flip_alleles_if_necessary ;
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
