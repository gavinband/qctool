
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_RELATEDNESS_COMPONENT_LAPACK_EIGEN_SOLVER
#define QCTOOL_RELATEDNESS_COMPONENT_LAPACK_EIGEN_SOLVER

#include <Eigen/Core>

namespace lapack
{
	void compute_eigendecomposition( Eigen::MatrixXd const& matrix, Eigen::VectorXd* eigenvalues, Eigen::MatrixXd* eigenvectors ) ;

	void compute_partial_eigendecomposition(
		Eigen::MatrixXd const& matrix,
		Eigen::VectorXd* eigenvalues,
		Eigen::MatrixXd* eigenvectors,
		double minimum_eigenvalue,
		double maximum_eigenvalue = std::numeric_limits< double >::max()
	) ;

	void compute_partial_eigendecomposition(
		Eigen::MatrixXd const& matrix,
		Eigen::VectorXd* eigenvalues,
		Eigen::MatrixXd* eigenvectors,
		int number_of_eigenvalues
	) ;
}

#endif
