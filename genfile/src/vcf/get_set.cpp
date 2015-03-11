
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set.hpp"

namespace genfile {
	namespace vcf {
		GenotypeSetterBase::GenotypeSetterBase( std::string const& scale ):
			m_probability_scale( scale == "identity" ? eIdentityScale : ePhredScale ),
			m_number_of_samples(0),
			m_sample(0),
			m_order_type( eUnorderedList ),
			m_value_type( eUnknownValueType ),
			m_number_of_entries(0),
			m_entry_i(0),
			m_missing( false )
		{
			if( scale != "identity" && scale != "phred" ) {
				throw BadArgumentError(
					"genfile::vcf::GenotypeSetterBase()",
					"scale=\"" + scale + "\"",
					"scale must be \"identity\" or \"phred\""
				) ;
			}
		}
		
		GenotypeSetterBase::~GenotypeSetterBase() throw() {}

		void GenotypeSetterBase::set_number_of_samples( std::size_t n ) {
			// destination is supposed to know its size.
			m_number_of_samples = n ;
		}

		void GenotypeSetterBase::set_number_of_alleles( std::size_t n ) {
			// handle only biallelic variants for now
			assert( n == 2 ) ;
		}

		bool GenotypeSetterBase::set_sample( std::size_t n ) {
			assert( n < m_number_of_samples ) ;
			// First test to see if no data was supplied (i.e. set_number_of_entries was uncalled) for the previous sample.
			//if( m_missing ) {
			//	set( n-1, 0.0, 0.0, 0.0 ) ;
			//}
			m_sample = n ;
			m_missing = true ;
			return true ;
		}

		void GenotypeSetterBase::set_order_type( OrderType const order_type, ValueType const value_type ) {
			assert(
				((value_type == eProbability) && (order_type == ePerUnorderedGenotype))
				|| 
				((value_type == eAlleleIndex) && (order_type == ePerOrderedHaplotype || order_type == ePerUnorderedHaplotype))
				||
				((value_type == eDosage) && (order_type == eBAlleleDosage))
			) ;
			m_order_type = order_type ;
			m_value_type = value_type ;
		}

		void GenotypeSetterBase::set_number_of_entries( std::size_t n ) {
			assert( n == 1 || n == 2 || n == 3 ) ;
			m_number_of_entries = n ;
			m_entry_i = 0 ;
			m_A = 0 ;
			m_B = 0 ;
			m_missing = false ;
		}

		void GenotypeSetterBase::operator()( MissingValue const value ) {
			// if any prob is missing, all are.
			m_missing = true ;
			if( ++m_entry_i == m_number_of_entries ) {
				set() ;
			}
		}

		void GenotypeSetterBase::set() {
			if( m_missing ) {
				set( m_sample, 0.0, 0.0, 0.0 ) ;
			} else if(
				(m_value_type == eDosage && m_order_type == eBAlleleDosage)
					|| ( m_value_type == eAlleleIndex && ( m_order_type == ePerOrderedHaplotype || m_order_type == ePerUnorderedHaplotype ))
			) {
				if( m_A == 0 && m_B == 2 ) {
					set( m_sample, 0.0, 0.0, 1.0 ) ;
				}
				else if( m_A == 1 && m_B == 1 ) {
					set( m_sample, 0.0, 1.0, 0.0 ) ;
				}
				else {
					set( m_sample, 1.0, 0.0, 0.0 ) ;
				}
			}
			else {
				set( m_sample, m_store[0], m_store[1], m_store[2] ) ;
			}
		}
		
		void GenotypeSetterBase::operator()( Integer const value ) {
			store( value ) ;
		}

		void GenotypeSetterBase::operator()( double const value ) {
			store( value ) ;
		}

		GenotypeSetter< SingleSNPGenotypeProbabilities >::GenotypeSetter( SingleSNPGenotypeProbabilities& result ):
			GenotypeSetterBase( "identity" ),
			m_result( result )
		{}

		void GenotypeSetter< SingleSNPGenotypeProbabilities >::set_number_of_samples( std::size_t n ) {
			// destination is supposed to know its size.
			m_result.resize( n ) ;
			GenotypeSetterBase::set_number_of_samples( n ) ;
		}

		void GenotypeSetter< SingleSNPGenotypeProbabilities >::set( std::size_t sample_i, double AA, double AB, double BB ) {
			m_result.set( sample_i, AA, AB, BB ) ;
		}

		GenotypeSetter< std::vector< double > >::GenotypeSetter( std::vector< double >& result ):
			GenotypeSetterBase( "identity" ),
			m_result( result )
		{}

		void GenotypeSetter< std::vector< double > >::set_number_of_samples( std::size_t n ) {
			m_result.resize( n * 3 ) ;
			GenotypeSetterBase::set_number_of_samples( n ) ;
		}

		void GenotypeSetter< std::vector< double > >::set( std::size_t sample_i, double AA, double AB, double BB ) {
			m_result[ sample_i * 3 + 0 ] = AA ;
			m_result[ sample_i * 3 + 1 ] = AB ;
			m_result[ sample_i * 3 + 2 ] = BB ;
		}
		
		ThreshholdingGenotypeSetter< std::vector< VariantEntry > >::
			ThreshholdingGenotypeSetter( std::vector< VariantEntry >& result, double threshhold ):
			GenotypeSetterBase( "identity" ),
			m_result( result ),
			m_threshhold( threshhold )
		{}

		void ThreshholdingGenotypeSetter< std::vector< VariantEntry > >::set_number_of_samples( std::size_t n ) {
			m_result.clear() ;
			m_result.resize( n, MissingValue() ) ;
			GenotypeSetterBase::set_number_of_samples( n ) ;
		}

		void ThreshholdingGenotypeSetter< std::vector< VariantEntry > >::set( std::size_t sample_i, double AA, double AB, double BB ) {
			if( AA > m_threshhold ) {
				m_result[ sample_i ] = 0 ;
			}
			else if( AB > m_threshhold ) {
				m_result[ sample_i ] = 1 ;
			}
			else if( BB > m_threshhold ) {
				m_result[ sample_i ] = 2 ;
			}
			else {
				m_result[ sample_i ] = MissingValue() ;
			}
		}

#if 0
		ThreshholdingGenotypeSetter< std::vector< int > >::
			ThreshholdingGenotypeSetter( std::vector< int >& result, double threshhold, int missing_value, int AA_value, int AB_value, int BB_value ):
			GenotypeSetterBase( "identity" ),
			m_result( result ),
			m_missing_value( missing_value ),
			m_AA_value( AA_value ),
			m_AB_value( AB_value ),
			m_BB_value( BB_value ),
			m_threshhold( threshhold )
		{}

		void ThreshholdingGenotypeSetter< std::vector< int > >::set_number_of_samples( std::size_t n ) {
			m_result.clear() ;
			m_result.resize( n, -1 ) ;
			GenotypeSetterBase::set_number_of_samples( n ) ;
		}

		void ThreshholdingGenotypeSetter< std::vector< int > >::set( std::size_t sample_i, double AA, double AB, double BB ) {
			if( AA > m_threshhold ) {
				m_result[ sample_i ] = m_AA_value ;
			}
			else if( AB > m_threshhold ) {
				m_result[ sample_i ] = m_AB_value ;
			}
			else if( BB > m_threshhold ) {
				m_result[ sample_i ] = m_BB_value ;
			}
			else {
				m_result[ sample_i ] = m_missing_value ;
			}
		}
#endif		
	}
}

