#include <Eigen/Core>
#include "LapackEigenDecomposition.hpp"
#include "../config.hpp"
#if HAVE_CLAPACK
	#include "clapack.h"
#else
extern "C" {
	/*
	*  function dsyev_ (see http://www.netlib.org/clapack/what/double/dsyev.c)
	*
	* Arguments:
	*
	JOBZ	(input) CHARACTER*1	  
			= 'N':	Compute eigenvalues only;	
			= 'V':	Compute eigenvalues and eigenvectors.	

	UPLO	(input) CHARACTER*1	  
			= 'U':	Upper triangle of A is stored;	 
			= 'L':	Lower triangle of A is stored.	 

	N		(input) INTEGER	  
			The order of the matrix A.	N >= 0.	  

	A		(input/output) DOUBLE PRECISION array, dimension (LDA, N)	
			On entry, the symmetric matrix A.  If UPLO = 'U', the	
			leading N-by-N upper triangular part of A contains the	 
			upper triangular part of the matrix A.	If UPLO = 'L',	 
			the leading N-by-N lower triangular part of A contains	 
			the lower triangular part of the matrix A.	 
			On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
			orthonormal eigenvectors of the matrix A.	
			If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')   
			or the upper triangle (if UPLO='U') of A, including the	  
			diagonal, is destroyed.	  

	LDA		(input) INTEGER	  
			The leading dimension of the array A.  LDA >= max(1,N).	  

	W		(output) DOUBLE PRECISION array, dimension (N)	 
			If INFO = 0, the eigenvalues in ascending order.   

	WORK	(workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
			On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

	LWORK	(input) INTEGER	  
			The length of the array WORK.  LWORK >= max(1,3*N-1).	
			For optimal efficiency, LWORK >= (NB+2)*N,	 
			where NB is the blocksize for DSYTRD returned by ILAENV.   

			If LWORK = -1, then a workspace query is assumed; the routine	
			only calculates the optimal size of the WORK array, returns	  
			this value as the first entry of the WORK array, and no error	
			message related to LWORK is issued by XERBLA.	

	INFO	(output) INTEGER   
			= 0:  successful exit	
			< 0:  if INFO = -i, the i-th argument had an illegal value	 
			> 0:  if INFO = i, the algorithm failed to converge; i	 
				  off-diagonal elements of an intermediate tridiagonal	 
				  form did not converge to zero.
	*/
	int dsyev_(char *jobz, char *uplo, int *n, double *a,
		 int *lda, double *w, double *work, int *lwork, 
		int *info) ;

	/*
	int dsyevr_(char *jobz, char *range, char *uplo, integer *n, 
		doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
		il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
		doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
		integer *lwork, integer *iwork, integer *liwork, integer *info)	// 

	doublereal dlamch_( char *cmach ) ;
	*/
}
#endif

#include <vector>

namespace lapack
{
	void compute_eigendecomposition( Eigen::MatrixXd const& matrix, Eigen::VectorXd* eigenvalues, Eigen::MatrixXd* eigenvectors ) {
		int N = matrix.cols() ;
		assert( matrix.rows() == N ) ;
		*eigenvectors = matrix ;
		int LDA = eigenvectors->outerStride() ;
		int LWORK = -1 ;
		int info = 0 ;
		char JOBZ = 'V' ;
		char UPLO = 'L' ;
		double work_size ;
		dsyev_( &JOBZ, &UPLO, &N, 0, &LDA, 0, &work_size, &LWORK, &info ) ;
		LWORK = work_size + 32 ;
		std::vector< double > workspace( work_size ) ;
		dsyev_( &JOBZ, &UPLO, &N, eigenvectors->data(), &LDA, eigenvalues->data(), &workspace[0], &LWORK, &info ) ;
		if( info != 0 ) {
			std::cerr << "!! compute_eigendecomposition(): info = " << info << ".\n" ;
			if( info < 0 ) {
				assert( 0 ) ;
			}
		}
	}

/*	void compute_eigendecomposition(
		Eigen::MatrixXd const& matrix,
		Eigen::VectorXd* eigenvalues,
		Eigen::MatrixXd* eigenvectors,
		double minimum_eigenvalue
	) {
		int N = matrix.cols() ;
		assert( matrix.rows() == N ) ;
		int M_LDA = matrix->outerStride() ;
		int EV_LDA = eigenvectors->outerStride() ;
		int LWORK = -1 ;
		int info = 0 ;
		char RANGE = 'V' ;
		char JOBZ = 'V' ;
		char UPLO = 'L' ;
		double maximum_eigenvalue = std::numeric_limits< double >::max() ;
		double ABSTOL = dlamch_( "Safe minimum" ) ;
		double work_size ;
		int M ;
		dsyev_( &JOBZ, &UPLO, &N, 0, &LDA, 0, &work_size, &LWORK, &info ) ;
		LWORK = work_size + 32 ;
		std::vector< double > workspace( work_size ) ;

		dsyevr_(
			&JOBZ, &RANGE, &UPLO, &N,
			matrix->data(),
			&M_LDA,
			&minimum_eigenvalue, &maximum_eigenvalue,
			0, 0,
			&ABSTOL,
			&M,
			eigenvalues->data(),
			eigenvectors->data(),
			&EV_LDA
			&workspace[0],
			&LWORK,
			&info
		) ;
		if( info != 0 ) {
			std::cerr << "!! compute_eigendecomposition(): info = " << info << ".\n" ;
			if( info < 0 ) {
				assert( 0 ) ;
			}
		}
	}*/
}
