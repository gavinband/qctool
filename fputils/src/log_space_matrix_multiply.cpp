
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include "fputils/floating_point_utils.hpp"
#include "fputils/log_space_matrix_multiply.hpp"

namespace fputils {
	Matrix log_space_matrix_multiply( Matrix const& left, Matrix const& right ) {
		assert( left.size2() == right.size1() ) ;
		Matrix result( left.size1(), right.size2() ) ;
		for( std::size_t i = 0; i < result.size1(); ++i ) {
			for( std::size_t j = 0; j < result.size2(); ++j ) {
				std::vector< double > entries( left.size2() ) ;
				for( std::size_t k = 0; k < left.size2(); ++k ) {
					entries[k] = left(i, k) + right( k, j ) ;
				}
				result( i , j ) = log_sum_exp( entries.begin(), entries.end() ) ;
			}
		}
		return result ;
	}
}
