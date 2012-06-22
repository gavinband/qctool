
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNP_IDENTIFYING_DATA2_HPP
#define GENFILE_SNP_IDENTIFYING_DATA2_HPP

#include <string>
#include <vector>
#include "genfile/GenomePosition.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	struct SNPIdentifyingData2
	{
	public:
		SNPIdentifyingData2() ;
		SNPIdentifyingData2(
			std::string const& SNPID,
			std::string const& RSID,
			GenomePosition const& position,
			char first_allele,
			char second_allele
		) ;
		SNPIdentifyingData2(
			std::string const& SNPID,
			std::string const& RSID,
			GenomePosition const& position,
			std::string const& first_allele,
			std::string const& second_allele
		) ;
		SNPIdentifyingData2(
			SNPIdentifyingData const& snp
		) ;
		SNPIdentifyingData2( SNPIdentifyingData2 const& other ) ;
		SNPIdentifyingData2& operator=( SNPIdentifyingData const& other ) ;
		SNPIdentifyingData2& operator=( SNPIdentifyingData2 const& other ) ;

		operator SNPIdentifyingData() const ;


		typedef string_utils::slice slice ;

		void set_SNPID( slice const& SNPID ) ;
		void set_rsid( slice const& rsid ) ;
		void set_position( GenomePosition const& position ) { m_position = position ;}
		void set_first_allele( slice const& allele ) ;
		void set_second_allele( slice const& allele ) ; 

		slice get_SNPID() const { return slice( m_data, 0, m_rsid_start ) ; }
		slice get_rsid() const { return slice( m_data, m_rsid_start, m_first_allele_start ) ; }
		GenomePosition const& get_position() const { return m_position ; }
		slice get_first_allele() const { return slice( m_data, m_first_allele_start, m_second_allele_start ) ; }
		slice get_second_allele() const { return slice( m_data, m_second_allele_start, m_data.size() ) ; }

		std::size_t get_estimated_bytes_used() const ;
	public:
		struct CompareFields {
			CompareFields( std::string const& fields_to_compare ) ;
			CompareFields( CompareFields const& other ) ;

			enum { eSNPID = 0x1, eRSID = 0x2, ePosition = 0x4, eAlleles = 0x8, eMask = 0xF } ;
			bool operator()( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) const ;
			bool are_equal( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) const ;
			bool check_if_comparable_fields_are_known( SNPIdentifyingData2 const& value ) const ;
		private:
			static std::vector< int > parse_fields_to_compare( std::string const& field_spec ) ;
			std::vector< int > const m_fields_to_compare ;
		} ;
	private:
		// we store all IDs and alleles in a string.
		// and hold slices to the string.
		std::string m_data ;
		typedef uint32_t Size ;
		// uint32_t m_SNPID_start ; // 0
		Size m_rsid_start ;
		Size m_first_allele_start ;
		Size m_second_allele_start ;
		GenomePosition m_position ;
		
		CompareFields& operator=( CompareFields const& other ) ;
	} ;	
	
	std::ostream& operator<<( std::ostream&, SNPIdentifyingData2 const& ) ;
	std::ostream& operator<<( std::ostream&, std::vector< SNPIdentifyingData2 > const& ) ;
	bool operator==( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) ;
	bool operator!=( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) ;
	bool operator<( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) ;
}

#endif
