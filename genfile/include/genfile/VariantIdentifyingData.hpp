
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNP_IDENTIFYING_DATA2_HPP
#define GENFILE_SNP_IDENTIFYING_DATA2_HPP

#include <string>
#include <vector>
#include <boost/function.hpp>
#include "genfile/GenomePosition.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	struct VariantIdentifyingData
	{
	public:
		VariantIdentifyingData( std::string const& rsid = "(unknown variant)" ) ;
		VariantIdentifyingData(
			std::string const& RSID,
			GenomePosition const& position,
			std::string const& first_allele,
			std::string const& second_allele
		) ;
	 	VariantIdentifyingData(
			std::string const& SNPID,
			std::string const& rsid,
			GenomePosition const& position,
			std::string const& first_allele,
			std::string const& second_allele
		) ;
		VariantIdentifyingData(
			SNPIdentifyingData const& snp
		) ;
		VariantIdentifyingData( VariantIdentifyingData const& other ) ;
		VariantIdentifyingData& operator=( VariantIdentifyingData const& other ) ;
		VariantIdentifyingData& operator=( SNPIdentifyingData const& other ) ;

		operator SNPIdentifyingData() const ;

		typedef string_utils::slice slice ;

		void clear_identifiers() ;
		std::size_t number_of_identifiers() const ;
		void set_primary_id( slice const& id ) ;
		void add_identifier( slice const& id ) ;
		std::vector< slice > get_identifiers(
			std::size_t start = 0,
			std::size_t end = std::string::npos
		) const ;
		void get_identifiers(
			boost::function< void( slice ) >,
			std::size_t start = 0,
			std::size_t end = std::string::npos
		) const ;
		std::string get_identifiers_as_string(
			std::string const& separator,
			std::size_t start = 0,
			std::size_t end = std::string::npos
		) const ;
		slice get_rsid() const { return slice( m_data, m_rsid_start, m_allele_starts[0] ) ; }

		void set_position( GenomePosition const& position ) { m_position = position ;}
		GenomePosition const& get_position() const { return m_position ; }

		void clear_alleles() ;
		void set_allele( std::size_t i, std::string const& allele ) ;
		void set_allele( std::size_t i, slice const& allele ) ;
		void set_allele( std::size_t i, char const* allele ) ;
		void add_allele( slice const& allele ) ;
		std::size_t number_of_alleles() const { return m_allele_starts.size() - 1 ; }
		slice get_allele( std::size_t i ) const ;
		std::vector< slice > get_alleles( std::size_t start = 0, std::size_t end = std::string::npos ) const ;
		void get_alleles( boost::function< void( slice ) >, std::size_t start = 0, std::size_t end = std::string::npos ) const ;
		void swap_alleles() ;

		std::size_t estimate_bytes_used() const ;
	public:
		struct CompareFields {
			CompareFields(
				std::string const& fields_to_compare = "position,rsid,IDs,alleles",
				bool flip_alleleles_if_necessary = false
			) ;
			CompareFields( CompareFields const& other ) ;
			CompareFields& operator=( CompareFields const& other ) ;

			bool get_flip_alleles_if_necessary() const { return m_flip_alleles_if_necessary ; }

			enum { eIDs = 0x1, eRSID = 0x2, ePosition = 0x4, eAlleles = 0x8, eMask = 0xF } ;
			bool operator()( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) const ;
			bool are_equal( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) const ;
			bool check_if_comparable_fields_are_known( VariantIdentifyingData const& value ) const ;
			std::vector< int > const& get_compared_fields() const { return m_fields_to_compare ;}
			
			std::string get_summary() const ;
			
		private:
			static std::vector< int > parse_fields_to_compare( std::string const& field_spec ) ;
			std::vector< int > m_fields_to_compare ;
			bool m_flip_alleles_if_necessary ;
		} ;
		
		friend bool operator==( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) ;
		friend bool operator!=( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) ;
		friend bool operator<( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) ;
		
	private:
		// we store all IDs and alleles in a string.
		// and hold slices to the string.
		std::string m_data ;
		typedef uint32_t Size ;
		Size m_rsid_start ; // always equals 0 in current implementation.
		std::vector< Size > m_allele_starts ;
		GenomePosition m_position ;
		
		CompareFields& operator=( CompareFields const& other ) ;
	} ;	
	
	std::ostream& operator<<( std::ostream&, VariantIdentifyingData const& ) ;
	std::ostream& operator<<( std::ostream&, std::vector< VariantIdentifyingData > const& ) ;
}

#endif
