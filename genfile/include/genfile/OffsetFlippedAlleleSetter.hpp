
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_FLIPPED_ALLELE_SETTER_HPP
#define GENFILE_FLIPPED_ALLELE_SETTER_HPP


#include <string>
#include <vector>
#include <algorithm>
#include <boost/format.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	struct OffsetFlippedAlleleSetter: public VariantDataReader::PerSampleSetter {
		static char const eUnknownFlip = '?' ;
		static char const eNoFlip = '+' ;
		static char const eFlip  = '-' ;

		~OffsetFlippedAlleleSetter() throw() {}
	
		OffsetFlippedAlleleSetter(
			VariantDataReader::PerSampleSetter& setter,
			std::size_t number_of_samples,
			char flip,
			std::size_t sample_offset
		):
			m_setter( setter ),
			m_number_of_samples( number_of_samples ),
			m_flip( flip ),
			m_sample_offset( sample_offset ),
			m_have_initialised( false ),
			m_number_of_alleles( 0 ),
			m_values( 3 ),
			m_entry_i( 0 )
		{}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			if( !m_have_initialised ) {
				m_number_of_alleles = nAlleles ;
				m_setter.initialise( m_number_of_samples, nAlleles ) ;
				m_have_initialised = true ;
			}
		}

		void set_offset( std::size_t offset ) {
			m_sample_offset = offset ;
		}

		std::size_t get_offset() const {
			return m_sample_offset ;
		}

		void set_flip( char flip ) {
            assert( flip == '+' || flip == '?' || flip == '-' ) ;
            m_flip = flip ;
        }
		std::size_t get_flip() const { return m_flip ; }
		
		bool set_sample( std::size_t n ) {
			assert( ( n + m_sample_offset ) < m_number_of_samples ) ;
			m_setter.set_sample( n + m_sample_offset ) ;
			return true ;
		}
		void set_number_of_entries( std::size_t n, OrderType order_type, ValueType value_type ) {
			m_values.resize( n ) ;
			m_order_type = order_type ;
			m_value_type = value_type ;
			m_entry_i = 0 ;
			m_setter.set_number_of_entries( n, order_type, value_type ) ;
		}

		void set_value( MissingValue const value ) { store( value ) ; }
		void set_value( std::string& value ) { store( value ) ; }
		void set_value( Integer const value ) { store( value ) ; }
		void set_value( double const value ) { store( value ) ; }

	private:
		VariantDataReader::PerSampleSetter& m_setter ;
		std::size_t const m_number_of_samples ;
		bool m_have_initialised ;
		char m_flip ;
		std::size_t m_sample_offset ;
		std::size_t m_number_of_alleles ;
		OrderType m_order_type ;
		ValueType m_value_type ;
		std::vector< VariantEntry > m_values ;
		std::size_t m_entry_i ;

		template< typename T >
		void store( T value ) {
			m_values[ m_entry_i++ ] = value ;
			if( m_entry_i == m_values.size() ) {
				set_values() ;
			}
		}

		void set_values() {
			if(
				m_flip == '+'
				|| !( m_order_type == ePerUnorderedGenotype || m_order_type == ePerAllele || m_value_type == eAlleleIndex )
			) {
				for( std::size_t i = 0; i < m_values.size(); ++i ) {
					m_setter.set_value( m_values[i] ) ;
				}
			} else if( m_flip == '-' ) {
				for( std::size_t i = 0; i < m_values.size(); ++i ) {
					VariantEntry entry ;
					std::size_t index = i ;
					if(
						m_order_type == ePerUnorderedGenotype
						|| m_order_type == ePerAllele
					) {
						// Values come one per allele or one per possible genotype.
						// Either way the right thing to do is reverse their order.
						index = m_values.size() - 1 - i ;
					}
					else if( m_order_type == ePerPhasedHaplotypePerAllele ) {
						// values come in the order p00 p01 ... p0(K-1) p10 p11 ... p1(K-1) ... p(J-1)(K-1)
						// where j=0,,(J-1) are the haplotypes and k=0..(K-1) are the alleles
						// Thus we can think of i=K*j+k and we want to output the value for index K*j+(K-1-k)
						index = ( m_number_of_alleles * ( i / m_number_of_alleles )) + ( m_number_of_alleles - 1 - ( i % m_number_of_alleles )) ;
					}
					else {
						// index = i
					}
				
					entry = m_values[ index ] ;
				
					// If values are allele indices they must be renumbered.
					if( !entry.is_missing() && m_value_type == eAlleleIndex ) {
						entry = VariantEntry::Integer( m_number_of_alleles - 1 ) - entry.as< VariantEntry::Integer >() ;
					}

					m_setter.set_value( entry ) ;
				}
			} else {
				// unknown flip.
				for( std::size_t i = 0; i < m_values.size(); ++i ) {
					m_setter.set_value( genfile::MissingValue() ) ;
				}
			}
		}
	} ;
}

#endif
