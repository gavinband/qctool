
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_INTERSECT_RANGES_HPP
#define CORE_INTERSECT_RANGES_HPP

#include <vector>
#include "metro/SampleRange.hpp"

namespace metro {
	namespace impl {
		// intersect two lists of ranges.
		// each range is assumed to be a sorted list of non-intersecting ranges of integers.
		// The result is also a sorted list of non-intersecting ranges.
		std::vector< metro::SampleRange > intersect_ranges(
			std::vector< metro::SampleRange > const& left,
			std::vector< metro::SampleRange > const& right
		) ;
	}
}

#endif
