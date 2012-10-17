#ifndef COMPUTE_X_XT_BLOCK_TASK_HPP
#define COMPUTE_X_XT_BLOCK_TASK_HPP

#include "../config.hpp"
#if HAVE_CBLAS
	#include "cblas.h"
#endif
#include <Eigen/Core>

//#define DEBUG_COMPUTE_XXT_BLOCK_TASK 1

namespace pca {
	struct ComputeXXtSymmetricBlockTask: public worker::Task {
		typedef Eigen::MatrixXd Matrix ;
		typedef Eigen::Block< Matrix > MatrixBlock ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::VectorBlock< Vector > VectorBlock ;

		ComputeXXtSymmetricBlockTask(
			MatrixBlock result,
			MatrixBlock non_missingness,
			VectorBlock const& data,
			VectorBlock const& non_missing_data,
			double scale
		):
			m_result( result ),
			m_non_missingness( non_missingness ),
			m_data( data ),
			m_non_missing_data( non_missing_data ),
			m_scale( scale )
		{}

	protected:
		MatrixBlock m_result ;
		MatrixBlock m_non_missingness ;
		Vector const m_data ;
		Vector const m_non_missing_data ;
		double const m_scale ;
	} ;

#if HAVE_CBLAS
	struct ComputeXXtSymmetricBlockUsingCblasTask: public ComputeXXtSymmetricBlockTask {
		ComputeXXtSymmetricBlockUsingCblasTask(
			MatrixBlock result,
			MatrixBlock non_missingness,
			VectorBlock const& data,
			VectorBlock const& non_missing_data,
			double scale
		):
			ComputeXXtSymmetricBlockTask( result, non_missingness, data, non_missing_data, scale )
		{}
		
		void operator()() {
			int const N = m_data.size() ;
			cblas_dsyr(
				CblasColMajor,
				CblasLower,
				N,
				m_scale,
				m_data.data(),
				m_data.innerStride(),
				m_result.data(),
				m_result.outerStride()
			) ;

			cblas_dsyr(
				CblasColMajor,
				CblasLower,
				N,
				1.0,
				m_non_missing_data.data(),
				m_non_missing_data.innerStride(),
				m_non_missingness.data(),
				m_non_missingness.outerStride()
			) ;
		}
	} ;
#endif	
	
	struct ComputeXXtSymmetricBlockUsingEigenTask: public ComputeXXtSymmetricBlockTask {
		ComputeXXtSymmetricBlockUsingEigenTask(
			MatrixBlock result,
			MatrixBlock non_missingness,
			VectorBlock const& data,
			VectorBlock const& non_missing_data,
			double scale
		):
			ComputeXXtSymmetricBlockTask( result, non_missingness, data, non_missing_data, scale )
		{}
		
		void operator()() {
			m_result
				.selfadjointView< Eigen::Lower >()
				.rankUpdate( m_data, m_scale ) ;
			m_non_missingness
				.selfadjointView< Eigen::Lower >()
				.rankUpdate( m_non_missing_data, 1.0 ) ;
		}
	} ;
	
	struct ComputeXXtBlockTask: public worker::Task {
		typedef Eigen::MatrixXd Matrix ;
		typedef Eigen::Block< Matrix > MatrixBlock ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::VectorBlock< Vector > VectorBlock ;

		ComputeXXtBlockTask(
			MatrixBlock result,
			MatrixBlock non_missingness,
			VectorBlock const& data1,
			VectorBlock const& non_missing_data1,
			VectorBlock const& data2,
			VectorBlock const& non_missing_data2,
			double scale
		):
			m_result( result ),
			m_non_missingness( non_missingness ),
			m_data1( data1 ),
			m_non_missing_data1( non_missing_data1 ),
			m_data2( data2 ),
			m_non_missing_data2( non_missing_data2 ),
			m_scale( scale )
		{}

	protected:
		MatrixBlock m_result ;
		MatrixBlock m_non_missingness ;
		Vector const m_data1 ;
		Vector const m_non_missing_data1 ;
		Vector const m_data2 ;
		Vector const m_non_missing_data2 ;
		double const m_scale ;
	} ;
	
#if HAVE_CBLAS
	struct ComputeXXtBlockUsingCblasTask: public ComputeXXtBlockTask {
		ComputeXXtBlockUsingCblasTask(
			MatrixBlock result,
			MatrixBlock non_missingness,
			VectorBlock const& data1,
			VectorBlock const& non_missing_data1,
			VectorBlock const& data2,
			VectorBlock const& non_missing_data2,
			double scale
		):
			ComputeXXtBlockTask( result, non_missingness, data1, non_missing_data1, data2, non_missing_data2, scale )
		{}
		
		void operator()() {
#if DEBUG_COMPUTE_XXT_BLOCK_TASK
			std::cerr << "m_data1.size() == " << m_data1.size() << ".\n" ;
			std::cerr << "m_data1.innerStride() == " << m_data1.innerStride() << ".\n" ;
			std::cerr << "m_data1.outerStride() == " << m_data1.outerStride() << ".\n" ;
			std::cerr << "m_data2.size() == " << m_data2.size() << ".\n" ;
			std::cerr << "m_result.innerStride() == " << m_result.innerStride() << ".\n" ;
			std::cerr << "m_result.outerStride() == " << m_result.outerStride() << ".\n" ;
#endif
			std::cerr.flush() ;

			{
				
				int const M = m_data1.size() ;
				int const N = m_data2.size() ;
/*
				int const lda = m_non_missing_data1.outerStride() ;
				int const ldb = m_non_missing_data2.outerStride() ;
*/				
				int const data1_stride = m_data1.innerStride() ;
				int const data2_stride = m_data2.innerStride() ;
				int const ldc = m_result.outerStride() ;
			
				cblas_dger(
					CblasColMajor,
					M, N,
					m_scale,
					m_data1.data(), data1_stride,
					m_data2.data(), data2_stride,
					m_result.data(),
					ldc
				) ;

/*
				cblas_dgemm(
					CblasColMajor,
					CblasNoTrans,
					CblasTrans,
					M, N, 1.0,
					m_scale,
					m_data1.data(), lda,
					m_data2.data(), ldb,
					1.0,
					m_result.data(),
					ldc
				) ;
*/
			}
			{
				int const M = m_non_missing_data1.size() ;
				int const N = m_non_missing_data2.size() ;
				int const data1_stride = m_non_missing_data1.innerStride() ;
				int const data2_stride = m_non_missing_data2.innerStride() ;
/*				
				int const lda = m_non_missing_data1.outerStride() ;
				int const ldb = m_non_missing_data2.outerStride() ;
*/
				int const ldc = m_non_missingness.outerStride() ;
				
				cblas_dger(
					CblasColMajor,
					M, N,
					m_scale,
					m_non_missing_data1.data(), data1_stride,
					m_non_missing_data2.data(), data2_stride,
					m_non_missingness.data(),
					ldc
				) ;
/*
				cblas_dgemm(
					CblasColMajor,
					CblasNoTrans,
					CblasTrans,
					M, N, 1.0,
					1.0,
					m_non_missing_data1.data(), lda,
					m_non_missing_data2.data(), ldb,
					1.0,
					m_non_missingness.data(),
					ldc
				) ;
*/
			}
		}
	} ;
#endif

	struct ComputeXXtBlockUsingEigenTask: public ComputeXXtBlockTask {
		ComputeXXtBlockUsingEigenTask(
			MatrixBlock result,
			MatrixBlock non_missingness,
			VectorBlock const& data1,
			VectorBlock const& non_missing_data1,
			VectorBlock const& data2,
			VectorBlock const& non_missing_data2,
			double scale
		):
			ComputeXXtBlockTask( result, non_missingness, data1, non_missing_data1, data2, non_missing_data2, scale )
		{}
		
		void operator()() {
			m_result += m_data1 * m_data2.transpose() * m_scale ;
			m_non_missingness += m_non_missing_data1 * m_non_missing_data2.transpose() ;
		}
	} ;
	
	
}

#endif
