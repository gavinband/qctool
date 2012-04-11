
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/vcf/get_set_eigen.hpp"

namespace genfile {
	namespace vcf {
		ThreshholdingGenotypeSetter< Eigen::VectorXd >::ThreshholdingGenotypeSetter(
			Eigen::VectorXd& result,
			double threshhold,
			double missing_value,
			double AA_value,
			double AB_value,
			double BB_value
		):
			m_result( result ),
			m_non_missingness( 0 ),
			m_missing_value( missing_value ),
			m_AA_value( AA_value ),
			m_AB_value( AB_value ),
			m_BB_value( BB_value ),
			m_threshhold( threshhold )
		{}

		ThreshholdingGenotypeSetter< Eigen::VectorXd >::ThreshholdingGenotypeSetter(
			Eigen::VectorXd& result,
			Eigen::VectorXd& non_missingness,
			double threshhold,
			double missing_value,
			double AA_value,
			double AB_value,
			double BB_value
		):
			m_result( result ),
			m_non_missingness( &non_missingness ),
			m_missing_value( missing_value ),
			m_AA_value( AA_value ),
			m_AB_value( AB_value ),
			m_BB_value( BB_value ),
			m_threshhold( threshhold )
		{}

		void ThreshholdingGenotypeSetter< Eigen::VectorXd >::set_number_of_samples( std::size_t n ) {
			m_result.setConstant( n, m_missing_value ) ;
			if( m_non_missingness ) {
				m_non_missingness->setConstant( n, 0 ) ;
			}
			GenotypeSetterBase::set_number_of_samples( n ) ;
		}

		void ThreshholdingGenotypeSetter< Eigen::VectorXd >::set( std::size_t sample_i, double AA, double AB, double BB ) {
			if( AA > m_threshhold ) {
				m_result( sample_i ) = m_AA_value ;
				if( m_non_missingness ) {
					(*m_non_missingness)( sample_i ) = 1 ;
				}
			}
			else if( AB > m_threshhold ) {
				m_result( sample_i ) = m_AB_value ;
				if( m_non_missingness ) {
					(*m_non_missingness)( sample_i ) = 1 ;
				}
			}
			else if( BB > m_threshhold ) {
				m_result( sample_i ) = m_BB_value ;
				if( m_non_missingness ) {
					(*m_non_missingness)( sample_i ) = 1 ;
				}
			}
			else {
				m_result( sample_i ) = m_missing_value ;
				if( m_non_missingness ) {
					(*m_non_missingness)( sample_i ) = 0 ;
				}
			}
		}
		
	}
}