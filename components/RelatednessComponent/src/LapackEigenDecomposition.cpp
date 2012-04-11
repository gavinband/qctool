
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <Eigen/Core>
#include "components/RelatednessComponent/LapackEigenDecomposition.hpp"
#include <iostream>
#include "../config.hpp"
#if HAVE_CLAPACK
	#include "clapack.h"
#else

extern "C" {
	/*
	*  function dsyev_ (see http://www.netlib.org/clapack/what/double/dsyev.c)

	-- LAPACK driver routine (version 3.0) --	
		Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,	 
		Courant Institute, Argonne National Lab, and Rice University	  
		June 30, 1999	


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
}

#if !HAVE_CLAPACK
extern "C" {
/*
 *	function dsyevr_ (see http://www.netlib.org/clapack/what/double/dsyevr.c)

	Arguments	
	=========	

	JOBZ	(input) CHARACTER*1	  
			= 'N':	Compute eigenvalues only;	
			= 'V':	Compute eigenvalues and eigenvectors.	

	RANGE	(input) CHARACTER*1	  
			= 'A': all eigenvalues will be found.	
			= 'V': all eigenvalues in the half-open interval (VL,VU]   
				   will be found.	
			= 'I': the IL-th through IU-th eigenvalues will be found.	
   ********* For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and	  
   ********* DSTEIN are called	 

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
			On exit, the lower triangle (if UPLO='L') or the upper	 
			triangle (if UPLO='U') of A, including the diagonal, is	  
			destroyed.	 

	LDA		(input) INTEGER	  
			The leading dimension of the array A.  LDA >= max(1,N).	  

	VL		(input) DOUBLE PRECISION   
	VU		(input) DOUBLE PRECISION   
			If RANGE='V', the lower and upper bounds of the interval to	  
			be searched for eigenvalues. VL < VU.	
			Not referenced if RANGE = 'A' or 'I'.	

	IL		(input) INTEGER	  
	IU		(input) INTEGER	  
			If RANGE='I', the indices (in ascending order) of the	
			smallest and largest eigenvalues to be returned.   
			1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.	
			Not referenced if RANGE = 'A' or 'V'.	

	ABSTOL	(input) DOUBLE PRECISION   
			The absolute error tolerance for the eigenvalues.	
			An approximate eigenvalue is accepted as converged	 
			when it is determined to lie in an interval [a,b]	
			of width less than or equal to	 

					ABSTOL + EPS *	 max( |a|,|b| ) ,	

			where EPS is the machine precision.	 If ABSTOL is less than	  
			or equal to zero, then	EPS*|T|	 will be used in its place,	  
			where |T| is the 1-norm of the tridiagonal matrix obtained	 
			by reducing A to tridiagonal form.	 

			See "Computing Small Singular Values of Bidiagonal Matrices	  
			with Guaranteed High Relative Accuracy," by Demmel and	 
			Kahan, LAPACK Working Note #3.	 

			If high relative accuracy is important, set ABSTOL to	
			DLAMCH( 'Safe minimum' ).  Doing so will guarantee that	  
			eigenvalues are computed to high relative accuracy when	  
			possible in future releases.  The current code does not	  
			make any guarantees about high relative accuracy, but	
			furutre releases will. See J. Barlow and J. Demmel,	  
			"Computing Accurate Eigensystems of Scaled Diagonally	
			Dominant Matrices", LAPACK Working Note #7, for a discussion   
			of which matrices define their eigenvalues to high relative	  
			accuracy.	

	M		(output) INTEGER   
			The total number of eigenvalues found.	0 <= M <= N.   
			If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.	  

	W		(output) DOUBLE PRECISION array, dimension (N)	 
			The first M elements contain the selected eigenvalues in   
			ascending order.   

	Z		(output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))	 
			If JOBZ = 'V', then if INFO = 0, the first M columns of Z	
			contain the orthonormal eigenvectors of the matrix A   
			corresponding to the selected eigenvalues, with the i-th   
			column of Z holding the eigenvector associated with W(i).	
			If JOBZ = 'N', then Z is not referenced.   
			Note: the user must ensure that at least max(1,M) columns are	
			supplied in the array Z; if RANGE = 'V', the exact value of M	
			is not known in advance and an upper bound must be used.   

	LDZ		(input) INTEGER	  
			The leading dimension of the array Z.  LDZ >= 1, and if	  
			JOBZ = 'V', LDZ >= max(1,N).   

	ISUPPZ	(output) INTEGER array, dimension ( 2*max(1,M) )   
			The support of the eigenvectors in Z, i.e., the indices	  
			indicating the nonzero elements in Z. The i-th eigenvector	 
			is nonzero only in elements ISUPPZ( 2*i-1 ) through	  
			ISUPPZ( 2*i ).	 
   ********* Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1   

	WORK	(workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
			On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

	LWORK	(input) INTEGER	  
			The dimension of the array WORK.  LWORK >= max(1,26*N).	  
			For optimal efficiency, LWORK >= (NB+6)*N,	 
			where NB is the max of the blocksize for DSYTRD and DORMTR	 
			returned by ILAENV.	  

			If LWORK = -1, then a workspace query is assumed; the routine	
			only calculates the optimal size of the WORK array, returns	  
			this value as the first entry of the WORK array, and no error	
			message related to LWORK is issued by XERBLA.	

	IWORK	(workspace/output) INTEGER array, dimension (LIWORK)   
			On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.	

	LIWORK	(input) INTEGER	  
			The dimension of the array IWORK.  LIWORK >= max(1,10*N).	

			If LIWORK = -1, then a workspace query is assumed; the	 
			routine only calculates the optimal size of the IWORK array,   
			returns this value as the first entry of the IWORK array, and	
			no error message related to LIWORK is issued by XERBLA.	  

	INFO	(output) INTEGER   
			= 0:  successful exit	
			< 0:  if INFO = -i, the i-th argument had an illegal value	 
			> 0:  Internal error   
*/
	int dsyevr_(
		char *jobz, char *range, char *uplo, int *n, 
		double *a, int *lda,
		double *vl, double *vu,
		int *il,
		int *iu,
		double *abstol,
		int *m, double *w, 
		double *z__, int *ldz,
		int *isuppz,
		double *work, int *lwork, int *iwork, int *liwork, int *info
	) ;

	double dlamch_( char *cmach ) ;
}
#endif

namespace lapack {
	void compute_eigendecomposition(
		Eigen::MatrixXd const& input_matrix,
		Eigen::VectorXd* eigenvalues,
		Eigen::MatrixXd* eigenvectors,
		double minimum_eigenvalue,
		double maximum_eigenvalue
	) {
		assert( input_matrix.rows() == input_matrix.cols() ) ;
		assert( eigenvalues ) ;
		assert( eigenvectors ) ;
		assert( minimum_eigenvalue == minimum_eigenvalue ) ;
		assert( maximum_eigenvalue == maximum_eigenvalue ) ;
		assert( minimum_eigenvalue < maximum_eigenvalue ) ;

		Eigen::MatrixXd matrix = input_matrix ;
		int N = matrix.cols() ;
		assert( matrix.rows() == N ) ;
		int M_LDA = matrix.outerStride() ;
		int EV_LDA = eigenvectors->outerStride() ;
		int info = 0 ;
		char RANGE = 'V' ;
		char JOBZ = 'V' ;
		char UPLO = 'L' ;
		double ABSTOL = dlamch_( const_cast< char* >( "Safe minimum" ) ) ;
		int number_of_eigenvalues ; // numbers of eigenvectors that are computed

		// set up workspaces
		int LWORK = -1 ;
		double work_size ;
		int LIWORK = -1 ;
		int iwork_size ;

		dsyevr_(
			&JOBZ, &RANGE, &UPLO, &N,
			matrix.data(), &M_LDA,
			&minimum_eigenvalue, &maximum_eigenvalue,
			0, 0,
			&ABSTOL,
			&number_of_eigenvalues, eigenvalues->data(),
			eigenvectors->data(), &EV_LDA,
			0,
			&work_size, &LWORK, &iwork_size, &LIWORK, &info
		) ;
		assert( info == 0 ) ;
		LWORK = work_size + 32 ;
		LIWORK = iwork_size + 32 ;
		std::vector< double > workspace( work_size ) ;
		std::vector< int > iworkspace( iwork_size ) ;

		// Now compute the decomposition.
		eigenvalues->resize( input_matrix.cols() ) ;
		eigenvectors->resize( input_matrix.cols(), input_matrix.cols() ) ;
		
		dsyevr_(
			&JOBZ, &RANGE, &UPLO, &N,
			const_cast< double* >( matrix.data() ), &M_LDA,
			&minimum_eigenvalue, &maximum_eigenvalue,
			0, 0,
			&ABSTOL,
			&number_of_eigenvalues, eigenvalues->data(),
			eigenvectors->data(), &EV_LDA,
			0,
			&workspace[0],
			&LWORK,
			&iworkspace[0],
			&LIWORK,
			&info
		) ;
		if( info != 0 ) {
			std::cerr << "!! compute_eigendecomposition(): info = " << info << ".\n" ;
			if( info < 0 ) {
				assert( 0 ) ;
			}
		} else {
			eigenvalues->resize( number_of_eigenvalues ) ;
			eigenvectors->resize( eigenvectors->rows(), number_of_eigenvalues ) ;
		}
	}
}

namespace lapack {
	void compute_eigendecomposition(
		Eigen::MatrixXd const& input_matrix,
		Eigen::VectorXd* eigenvalues,
		Eigen::MatrixXd* eigenvectors,
		int number_of_eigenvalues
	) {
		assert( input_matrix.rows() == input_matrix.cols() ) ;
		assert( number_of_eigenvalues > 0 ) ;
		assert( eigenvalues ) ;
		assert( eigenvectors ) ;
		assert( number_of_eigenvalues <= input_matrix.cols() ) ;
		
		Eigen::MatrixXd matrix = input_matrix ;
		int N = matrix.cols() ;
		assert( matrix.rows() == N ) ;
		int M_LDA = matrix.outerStride() ;
		int EV_LDA = eigenvectors->outerStride() ;
		int info = 0 ;
		char RANGE = 'I' ;
		char JOBZ = 'V' ;
		char UPLO = 'L' ;
		double ABSTOL = dlamch_( const_cast< char* >( "Safe minimum" ) ) ;
		int IL = 1;
		int IU = number_of_eigenvalues + 1 ;

		// set up workspaces
		int LWORK = -1 ;
		double work_size ;
		int LIWORK = -1 ;
		int iwork_size ;

		dsyevr_(
			&JOBZ, &RANGE, &UPLO, &N,
			matrix.data(), &M_LDA,
			0, 0,
			&IL, &IU,
			&ABSTOL,
			&number_of_eigenvalues, eigenvalues->data(),
			eigenvectors->data(), &EV_LDA,
			0,
			&work_size, &LWORK, &iwork_size, &LIWORK, &info
		) ;
		assert( info == 0 ) ;
		LWORK = work_size + 32 ;
		LIWORK = iwork_size + 32 ;
		std::vector< double > workspace( work_size ) ;
		std::vector< int > iworkspace( iwork_size ) ;

		// compute the decomposition
		eigenvalues->resize( number_of_eigenvalues ) ;
		eigenvectors->resize( input_matrix.cols(), number_of_eigenvalues ) ;
		
		dsyevr_(
			&JOBZ, &RANGE, &UPLO, &N,
			const_cast< double* >( matrix.data() ), &M_LDA,
			0, 0,
			&IL, &IU,
			&ABSTOL,
			&number_of_eigenvalues, eigenvalues->data(),
			eigenvectors->data(), &EV_LDA,
			0,
			&workspace[0],
			&LWORK,
			&iworkspace[0],
			&LIWORK,
			&info
		) ;
		if( info != 0 ) {
			std::cerr << "!! compute_eigendecomposition(): info = " << info << ".\n" ;
			if( info < 0 ) {
				assert( 0 ) ;
			}
		} else {
			eigenvalues->resize( number_of_eigenvalues ) ;
			eigenvectors->resize( eigenvectors->rows(), number_of_eigenvalues ) ;
		}
	}
}
