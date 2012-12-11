
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
			GenotypeSetter( Eigen::MatrixBase< Derived >& result ): m_result( result ) {} ;
			void set_number_of_samples( std::size_t n ) {
				m_result.derived().resize( n, 3 ) ;
				GenotypeSetterBase::set_number_of_samples( n ) ;
			}
			void set( std::size_t sample_i, double AA, double AB, double BB ) {
				m_result( sample_i, 0 ) = AA ;
				m_result( sample_i, 1 ) = AB ;
				m_result( sample_i, 2 ) = BB ;
			}
		private:
			Eigen::MatrixBase< Derived >& m_result ;
		} ;
	
		// Genotype setter which stores hard genotype calls as values in an Eigen::Vector with specified values for each genotype.
		template<>
		struct ThreshholdingGenotypeSetter< Eigen::VectorXd >: public GenotypeSetterBase
		{
			ThreshholdingGenotypeSetter(
				Eigen::VectorXd& result,
				double threshhold
			) ;
			ThreshholdingGenotypeSetter(
				Eigen::VectorXd& result,
				double threshhold,
				double missing_value,
				double AA_value,
				double AB_value,
				double BB_value
			) ;
			ThreshholdingGenotypeSetter(
				Eigen::VectorXd& result,
				Eigen::VectorXd& non_missingness,
				double threshhold,
				double missing_value = -1,
				double AA_value = 0,
				double AB_value = 1,
				double BB_value = 2
			) ;
			void set_number_of_samples( std::size_t n ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			Eigen::VectorXd& m_result ;
			Eigen::VectorXd* m_non_missingness ;
			int const m_missing_value ;
			int const m_AA_value ;
			int const m_AB_value ;
			int const m_BB_value ;
			double const m_threshhold ;
		} ;
		
		// Genotype setter which stores hard genotype calls as values in an Eigen::Vector with specified values for each genotype.
		template< typename HaplotypeAlleles >
		struct PhasedGenotypeSetter: public VariantDataReader::PerSampleSetter
		{
			typedef HaplotypeAlleles NonMissingness ;
			PhasedGenotypeSetter(
				HaplotypeAlleles& result,
				double missing_value = std::numeric_limits< double >::quiet_NaN()
			):
				m_result( result ),
				m_non_missingness( 0 ),
				m_missing_value( missing_value )
			{}

			PhasedGenotypeSetter(
				HaplotypeAlleles& result,
				NonMissingness& non_missingness,
				double missing_value = 0
			):
				m_result( result ),
				m_non_missingness( &non_missingness ),
				m_missing_value( missing_value )
			{}

			void set_number_of_samples( std::size_t n ) {
				m_result.resize( n, 2  ) ;
				if( m_non_missingness ) {
					m_non_missingness->resize( n, 2 ) ;
				}
				m_sample_i = 0 ;
				m_entry_i = 0 ;
			}

			void set_sample( std::size_t i ) {
				m_sample_i = i ;
				m_entry_i = 0 ;
			}

			void set_number_of_entries( std::size_t n ) {
				if( n != 2 ) {
					throw BadArgumentError( "genfile::vcf::PhasedGenotypeSetter::set_number_of_entries()", "Only a value of 2 is supported." ) ;
				}
			}

			void set_order_type( OrderType const type ) {
				if( type != eOrderedList ) {
					throw BadArgumentError( "genfile::vcf::PhasedGenotypeSetter::set_order_type()", "Only phased genotypes (an ordered list) is supported." ) ;
				}
			}

			void operator()( MissingValue const value ) {
				m_result( m_sample_i, m_entry_i ) = m_missing_value ;
				if( m_non_missingness ) {
					(*m_non_missingness)( m_sample_i, m_entry_i ) = 0 ;
				}
				++m_entry_i ;
			}

			void operator()( Integer const value ) {
				m_result( m_sample_i, m_entry_i ) = value ;
				if( m_non_missingness ) {
					(*m_non_missingness)( m_sample_i, m_entry_i ) = 1 ;
				}
				++m_entry_i ;
			}

		private:
			HaplotypeAlleles& m_result ;
			NonMissingness* m_non_missingness ;
			int const m_missing_value ;
			int m_sample_i ;
			int m_entry_i ;
		} ;
	}
}

#endif
