#ifndef SNPTEST_FINITELY_SUPPORTED_FUNCTION_SET_HPP
#define SNPTEST_FINITELY_SUPPORTED_FUNCTION_SET_HPP

#include <memory>
#include "Eigen/Core"

namespace snptest {
	// struct FinitelySupportedFunctionSet
	// FinitelySupportedFunctionSet is a value type that represents a collection
	// of functions, all having support contained within the same finite set
	// of points on the real line.
	// An example would be a numerical-valued covariate whose values for each sample
	// are only known with some uncertainty (so we only have a distribution.) 
	struct FinitelySupportedFunctionSet {
	public:
		typedef std::auto_ptr< FinitelySupportedFunctionSet > UniquePtr ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
		typedef Matrix::RowXpr RowExpression ;
		typedef Matrix::ConstRowXpr ConstRowExpression ;
		typedef Matrix::ColXpr ColumnExpression ;
		typedef Matrix::ConstColXpr ConstColumnExpression ;
		
		FinitelySupportedFunctionSet(
			Vector const& support,
			Matrix const& matrix
		) ;

		FinitelySupportedFunctionSet(
			FinitelySupportedFunctionSet const& other
		) ;

		FinitelySupportedFunctionSet& operator=(
			FinitelySupportedFunctionSet const& other
		) ;
		
		std::size_t get_size_of_support() const {
			return m_support.size() ;
		}

		Vector const& get_support() const {
			return m_support ;
		}

		double const& get_support( std::size_t i ) const {
			return m_support(i) ;
		}

		double& get_support( std::size_t i ) {
			return m_support(i) ;
		}

		std::size_t get_number_of_functions() const {
			return m_matrix.rows() ;
		}
		
		ConstRowExpression get_values( std::size_t i ) const {
			return m_matrix.row( i ) ;
		}

		RowExpression get_values( std::size_t i ) {
			return m_matrix.row( i ) ;
		}
		
		Matrix const& get_values() const {
			return m_matrix ;
		}
		
	protected:
		Vector m_support ;
		Matrix m_matrix ;
	} ;
}

#endif
