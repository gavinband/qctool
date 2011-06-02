#include <memory>
#include <cassert>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "snptest/FinitelySupportedFunctionSet.hpp"

namespace snptest {
	FinitelySupportedFunctionSet::FinitelySupportedFunctionSet(
		FinitelySupportedFunctionSet::Vector const& support,
		Matrix const& matrix
	):
		m_support( support ),
		m_matrix( matrix )
	{
		assert( m_support.size() == m_matrix.cols() ) ;
	}
}
