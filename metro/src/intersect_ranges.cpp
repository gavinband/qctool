
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include "metro/SampleRange.hpp"
#include "metro/intersect_ranges.hpp"

namespace metro {
	namespace impl {
		std::vector< metro::SampleRange > intersect_ranges(
			std::vector< metro::SampleRange > const& left,
			std::vector< metro::SampleRange > const& right
		) {
			std::vector< metro::SampleRange > result ;
			std::size_t i = 0, j = 0 ;
			while( i < left.size() && j < right.size() ) {
				int range_start = std::max( left[i].begin(), right[j].begin() ) ;
				int range_end = std::min( left[i].end(), right[j].end() ) ;
				if( range_end > range_start ) {
					result.push_back( metro::SampleRange( range_start, range_end )) ;
				}
				if( left[i].end() == range_end ) {
					++i ;
				}
				if( right[j].end() == range_end ) {
					++j ;
				}
			}
			return result ;
		}
	}
}
