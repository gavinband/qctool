
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
		VariantIdentifyingData() ;
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

		void set_rsid( slice const& rsid ) ;
		void set_genome_position( GenomePosition const& position ) { m_position = position ;}
		void set_chromosome( Chromosome const& chromosome ) { m_position.chromosome() = chromosome ;}
		void set_first_allele( slice const& allele ) ;
		void set_second_allele( slice const& allele ) ;
		void add_allele( slice const& allele ) ;

		void clear_identifiers() ;
		void add_identifier( slice const& id ) ;
		void swap_alleles() ;
		
		slice get_rsid() const { return slice( m_data, m_rsid_start, m_allele_starts[0] ) ; }
		
		GenomePosition const& get_position() const { return m_position ; }
		std::size_t number_of_alleles() const { return m_allele_starts.size() ; }
		slice get_allele( std::size_t i ) const {
			assert( i < m_allele_starts.size() ) ;
			if( i+1 < m_allele_starts.size() ) {
				return slice( m_data, m_allele_starts[i], m_allele_starts[i+1] ) ;
			} else {
				return slice( m_data, m_allele_starts[i], m_identifiers_start ) ;
			}
		}
		void get_alleles( boost::function< void( slice ) > ) const ;

		std::vector< slice > get_alternative_identifiers() const ;
		void get_alternative_identifiers( boost::function< void( slice ) > ) const ;
		std::string get_alternate_identifiers_as_string() const ;

		std::size_t get_estimated_bytes_used() const ;
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
		Size m_identifiers_start ;
		GenomePosition m_position ;
		
		CompareFields& operator=( CompareFields const& other ) ;
	} ;	
	
	std::ostream& operator<<( std::ostream&, VariantIdentifyingData const& ) ;
	std::ostream& operator<<( std::ostream&, std::vector< VariantIdentifyingData > const& ) ;
}

#endif
