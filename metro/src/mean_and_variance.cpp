
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <Eigen/Dense>
#include "metro/mean_and_variance.hpp"

namespace metro {
	OnlineElementwiseMeanAndVariance::Storage OnlineElementwiseMeanAndVariance::get_mean() const {
		// Return m_mean, but with NaN for counts that are zero.
		return m_mean.array() + (( m_nonmissingness.array() / m_nonmissingness.array() ) - 1 ) ;
	}

	OnlineElementwiseMeanAndVariance::Storage OnlineElementwiseMeanAndVariance::get_variance() const {
		Storage result = m_sum_of_squares_of_differences.array() / ( m_nonmissingness.array() - 1 ) ;
		// Set result to NaN when there is 0 or 1 observation.
		result.array() *= m_nonmissingness.array() / m_nonmissingness.array() ;
		result.array() *= ( m_nonmissingness.array() - 1 ) / ( m_nonmissingness.array() - 1 ) ;
		return result ;
	}
	
	double OnlineElementwiseMeanAndVariance::get_count( int row, int column ) const {
		return m_nonmissingness(row, column) ;
	}

	double OnlineElementwiseMeanAndVariance::get_mean( int row, int column ) const {
		// Return m_mean, but with NaN for counts that are zero.
		return m_mean(row, column) + (( m_nonmissingness(row, column) / m_nonmissingness(row, column) ) - 1 ) ;
	}

	double OnlineElementwiseMeanAndVariance::get_variance( int row, int column ) const {
		// The result is NaN unless there are at least 2 observations
		double result = std::numeric_limits< double >::quiet_NaN() ;
		if( m_nonmissingness(row, column) > 1 ) {
			result = m_sum_of_squares_of_differences(row, column) / ( m_nonmissingness(row, column) - 1 ) ;
		}
		return result ;
	}
}
