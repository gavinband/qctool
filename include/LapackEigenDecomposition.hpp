#ifndef QCTOOL_LAPACK_EIGEN_SOLVER
#define QCTOOL_LAPACK_EIGEN_SOLVER

#include <Eigen/Core>

namespace lapack
{
	void compute_eigendecomposition( Eigen::MatrixXd const& matrix, Eigen::VectorXd* eigenvalues, Eigen::MatrixXd* eigenvectors ) ;
}

#endif
