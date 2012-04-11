
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_HAPLOTYPE_FREQUENCY_LOGLIKELIHOOD_HPP
#define QCTOOL_HAPLOTYPE_FREQUENCY_LOGLIKELIHOOD_HPP

#include <vector>
#include <Eigen/Core>

// Loglikelihood of table of genotypes at 2 SNPs, given parameters specifying the haplotype frequencies.
// The parameters are \pi_01, \pi_10, and \pi_11.  Then \pi_00 is one minus the sum of the others.
// \pi_ab is the frequency of the haplotype with genotype a at SNP 1 and b at SNP 2.
// Now the probability of any table is obtained
// 
struct HaplotypeFrequencyLogLikelihood {
	typedef Eigen::MatrixXd Matrix ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::RowVectorXd RowVector ;

	HaplotypeFrequencyLogLikelihood( Matrix const& genotype_table ) ;
	void evaluate_at( Vector const& pi ) ;
	double get_value_of_function() const ;
	Vector get_value_of_first_derivative() ;
	Matrix get_value_of_second_derivative() ;

	Vector get_MLE_by_EM() const ;

	private:
		Matrix const m_genotype_table ;
		std::vector< std::vector< RowVector > > m_dpi ;
		double m_ll ;
		RowVector m_D_ll ;
		Matrix m_DDt_ll ;
	private:
		Vector estimate_parameters( double const AB_ab ) const ;
} ;

#endif

