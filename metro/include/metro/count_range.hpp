
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_COUNT_RANGES_HPP
#define CORE_COUNT_RANGES_HPP

#include <vector>
#include "metro/SampleRange.hpp"

namespace metro {
	namespace impl {
		// count the number of samples in a range
		// each list is assumed to be a sorted list of non-intersecting ranges of integers.
		int count_range(
			std::vector< metro::SampleRange > const& range
		) ;
	}
}

#endif
