
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
	
	FinitelySupportedFunctionSet::FinitelySupportedFunctionSet(
		FinitelySupportedFunctionSet const& other
	):
		m_support( other.m_support ),
		m_matrix( other.m_matrix )
	{}

	FinitelySupportedFunctionSet& FinitelySupportedFunctionSet::operator=(
		FinitelySupportedFunctionSet const& other
	) {
		m_support = other.m_support ;
		m_matrix = other.m_matrix ;
		return *this ;
	}
	
}
