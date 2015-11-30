
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_GET_SET_EIGEN_HPP
#define GENFILE_VCF_GET_SET_EIGEN_HPP

#include <Eigen/Core>
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	namespace vcf {
		// Genotype setter which stores its genotypes in an Eigen matrix with
		// 3 rows and n columns (where n is the number of samples.)
		template< typename Derived >
		struct GenotypeSetter< Eigen::MatrixBase< Derived > >: public GenotypeSetterBase
		{
			GenotypeSetter( Eigen::MatrixBase< Derived >& result, std::string const& scale = "identity" ): GenotypeSetterBase( scale ), m_result( result ) {
				result.setZero() ;
			} ;
			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				assert( nAlleles ==  2 ) ;
				m_result.derived().resize( nSamples, 3 ) ;
				GenotypeSetterBase::initialise( nSamples, nAlleles ) ;
			}
			void set( std::size_t sample_i, double AA, double AB, double BB ) {
				m_result( sample_i, 0 ) = AA ;
				m_result( sample_i, 1 ) = AB ;
				m_result( sample_i, 2 ) = BB ;
			}
		private:
			Eigen::MatrixBase< Derived >& m_result ;
		} ;
	
		namespace impl {
			template< typename Data >
			struct ThreshholdedCallGetter1: public GenotypeSetterBase {
				ThreshholdedCallGetter1( ThreshholdedCallGetter1 const& other ):
					m_result( other.m_result ),
					m_threshhold( other.m_threshhold ),
					m_missing_value( other.m_missing_value ),
					m_AA_value( other.m_AA_value ),
					m_AB_value( other.m_AB_value ),
					m_BB_value( other.m_BB_value )
				{}

				ThreshholdedCallGetter1(
					Data& result,
					double threshhold,
					double missing_value = -1,
					double AA_value = 0,
					double AB_value = 1,
					double BB_value = 2
				):
					m_result( result ),
					m_threshhold( threshhold ),
					m_missing_value( missing_value ),
					m_AA_value( AA_value ),
					m_AB_value( AB_value ),
					m_BB_value( BB_value )
				{}

				void initialise( std::size_t nSamples, std::size_t nAlleles ) {
					//m_result.setConstant( n, m_missing_value ) ;
					assert( nAlleles == 2 ) ;
					m_result.resize( nSamples ) ;
					GenotypeSetterBase::initialise( nSamples, nAlleles ) ;
				}

				void set( std::size_t sample_i, double AA, double AB, double BB ) {
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
			private:
				Data& m_result ;
				double const m_threshhold ;
				int const m_missing_value ;
				int const m_AA_value ;
				int const m_AB_value ;
				int const m_BB_value ;
			} ;

			template< typename Data, typename Nonmissingness >
			struct ThreshholdedCallGetter2: public GenotypeSetterBase {
				ThreshholdedCallGetter2( ThreshholdedCallGetter2 const& other ):
					m_result( other.m_result ),
					m_nonmissingness( other.m_nonmissingness ),
					m_threshhold( other.m_threshhold ),
					m_missing_value( other.m_missing_value ),
					m_AA_value( other.m_AA_value ),
					m_AB_value( other.m_AB_value ),
					m_BB_value( other.m_BB_value )
				{}
				
				ThreshholdedCallGetter2(
					Data& result,
					Nonmissingness& nonmissingness,
					double threshhold,
					double missing_value = -1,
					double AA_value = 0,
					double AB_value = 1,
					double BB_value = 2
				):
					m_result( result ),
					m_nonmissingness( nonmissingness ),
					m_threshhold( threshhold ),
					m_missing_value( missing_value ),
					m_AA_value( AA_value ),
					m_AB_value( AB_value ),
					m_BB_value( BB_value )
				{}

				void initialise( std::size_t nSamples, std::size_t nAlleles ) {
					assert( nAlleles == 2 ) ;
					m_result.setConstant( nSamples, m_missing_value ) ;
					m_nonmissingness.setConstant( nSamples, 0 ) ;
					GenotypeSetterBase::initialise( nSamples, nAlleles ) ;
				}

				void set( std::size_t sample_i, double AA, double AB, double BB ) {
					if( AA > m_threshhold ) {
						m_result( sample_i ) = m_AA_value ;
						m_nonmissingness( sample_i ) = 1 ;
					}
					else if( AB > m_threshhold ) {
						m_result( sample_i ) = m_AB_value ;
						m_nonmissingness( sample_i ) = 1 ;
					}
					else if( BB > m_threshhold ) {
						m_result( sample_i ) = m_BB_value ;
						m_nonmissingness( sample_i ) = 1 ;
					}
					else {
						m_result( sample_i ) = m_missing_value ;
						m_nonmissingness( sample_i ) = 0 ;
					}
				}
			private:
				Data& m_result ;
				Nonmissingness& m_nonmissingness;
				double const m_threshhold ;
				int const m_missing_value ;
				int const m_AA_value ;
				int const m_AB_value ;
				int const m_BB_value ;
			} ;
		}

		template< typename Data >
		impl::ThreshholdedCallGetter1< Data > get_threshholded_calls(
			Data& result,
			double threshhold,
			double missing_value = -1, 
			double AA_value = 0,
			double AB_value = 1,
			double BB_value = 2
		) {
			return impl::ThreshholdedCallGetter1< Data >( result, threshhold, missing_value, AA_value, AB_value, BB_value ) ;
		}
 
		template< typename Data, typename Nonmissingness >
		impl::ThreshholdedCallGetter2< Data, Nonmissingness > get_threshholded_calls(
			Data& result,
			Nonmissingness& nonmissingness,
			double threshhold,
			double missing_value = -1, 
			double AA_value = 0,
			double AB_value = 1,
			double BB_value = 2
		) {
			return impl::ThreshholdedCallGetter2< Data, Nonmissingness >( result, nonmissingness, threshhold, missing_value, AA_value, AB_value, BB_value ) ;
		}

		// Genotype setter which stores hard genotype calls as values in an Eigen::Vector with specified values for each genotype.
		template< typename HaplotypeAlleles >
		struct PhasedGenotypeSetter: public VariantDataReader::PerSampleSetter
		{
			typedef HaplotypeAlleles NonMissingness ;
			~PhasedGenotypeSetter() throw() {}
			PhasedGenotypeSetter(
				HaplotypeAlleles& result,
				double missing_value = std::numeric_limits< double >::quiet_NaN()
			):
				m_result( result ),
				m_non_missingness( 0 ),
				m_missing_value( missing_value ),
				m_allele_coding( 2, 0 )
			{
				m_result.setConstant( m_missing_value ) ;
				m_allele_coding[0] = 0 ;
				m_allele_coding[1] = 1 ;
			}

			PhasedGenotypeSetter(
				HaplotypeAlleles& result,
				NonMissingness& non_missingness,
				double missing_value = 0,
				double const A_coding = 0,
				double const B_coding = 1
			):
				m_result( result ),
				m_non_missingness( &non_missingness ),
				m_missing_value( missing_value ),
				m_allele_coding( 2, 0 )
			{
				m_result.setConstant( m_missing_value ) ;
				m_allele_coding[0] = A_coding ;
				m_allele_coding[1] = B_coding ;
			}

			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				assert( nAlleles == 2 ) ;
				m_result.resize( nSamples, 2  ) ;
				if( m_non_missingness ) {
					m_non_missingness->resize( nSamples, 2 ) ;
				}
				m_sample_i = 0 ;
				m_entry_i = 0 ;
			}

			void set_number_of_alleles( std::size_t n ) {
				assert( n == 2 ) ;
			}

			bool set_sample( std::size_t i ) {
				m_sample_i = i ;
				m_entry_i = 0 ;
				return true ;
			}

			void set_number_of_entries( uint32_t, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				if( n != 2 ) {
					throw BadArgumentError( "genfile::vcf::PhasedGenotypeSetter::set_number_of_entries()", "Only a value of 2 is supported." ) ;
				}
				assert( order_type == ePerOrderedHaplotype ) ;
				assert( value_type == eAlleleIndex ) ;
			}

			void set_order_type( OrderType const type ) {
				if( type != eOrderedList ) {
					throw BadArgumentError( "genfile::vcf::PhasedGenotypeSetter::set_order_type()", "Only phased genotypes (an ordered list) is supported." ) ;
				}
			}

			void set_value( MissingValue const value ) {
				m_result( m_sample_i, m_entry_i ) = m_missing_value ;
				if( m_non_missingness ) {
					(*m_non_missingness)( m_sample_i, m_entry_i ) = 0 ;
				}
				++m_entry_i ;
			}

			void set_value( Integer const value ) {
				assert( value >= 0 && value < 2 ) ;
				m_result( m_sample_i, m_entry_i ) = m_allele_coding[ std::size_t( value ) ] ;
				if( m_non_missingness ) {
					(*m_non_missingness)( m_sample_i, m_entry_i ) = 1 ;
				}
				++m_entry_i ;
			}

		private:
			HaplotypeAlleles& m_result ;
			NonMissingness* m_non_missingness ;
			double const m_missing_value ;
			std::vector< double > m_allele_coding ;
			int m_sample_i ;
			int m_entry_i ;
		} ;
	}
}

#endif
