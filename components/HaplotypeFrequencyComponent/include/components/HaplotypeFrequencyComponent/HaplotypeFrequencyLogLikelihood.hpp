
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_HAPLOTYPE_FREQUENCY_LOGLIKELIHOOD_HPP
#define QCTOOL_HAPLOTYPE_FREQUENCY_LOGLIKELIHOOD_HPP

#include <vector>
#include <Eigen/Core>
#include <boost/math/distributions/binomial.hpp>

// Loglikelihood of table of genotypes  and known haplotypes at 2 biallelic variants,
// given parameters specifying the haplotype frequencies.
// The parameters are \pi_01, \pi_10, and \pi_11.  Then \pi_00 can be
// computed as one minus the sum of the others.
//
// This loglikelihood can additionally estimate the parameters using the EM algorithm.
//
struct HaplotypeFrequencyLogLikelihood {
	typedef Eigen::MatrixXd Matrix ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::RowVectorXd RowVector ;

	HaplotypeFrequencyLogLikelihood(
		Matrix const& genotype_table,
		Matrix const& haplotype_table = Matrix::Zero(2,2)
	) ;

	void evaluate_at( Vector const& pi ) ;
	double get_value_of_function() const ;
	Vector get_value_of_first_derivative() ;
	Matrix get_value_of_second_derivative() ;

	Vector get_MLE_by_EM() const ;

	private:
		Matrix const m_genotype_table ;
		Matrix const m_haplotype_table ;
		std::vector< std::vector< RowVector > > m_dpi ;
		double m_ll ;
		RowVector m_D_ll ;
		Matrix m_DDt_ll ;
	private:
		Vector estimate_parameters( double const p ) const ;
} ;

#endif

