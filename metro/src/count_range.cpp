
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <algorithm>
#include "metro/SampleRange.hpp"
#include "metro/count_range.hpp"

namespace metro {
	namespace impl {
		int count_range(
			std::vector< metro::SampleRange > const& range
		) {
			int result = 0 ;
			for( std::size_t i = 0; i < range.size(); ++i ) {
				result += range[i].size() ;
			}
			return result ;
		}
	}
}

