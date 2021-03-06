
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_UNION_RANGES_HPP
#define CORE_UNION_RANGES_HPP

#include <vector>
#include "metro/SampleRange.hpp"

namespace metro {
	namespace impl {
		// union two lists of ranges.
		// each list is assumed to be a sorted list of non-intersecting ranges of integers.
		// The result is also a sorted list of non-intersecting ranges.
		std::vector< metro::SampleRange > union_ranges(
			std::vector< metro::SampleRange > const& left,
			std::vector< metro::SampleRange > const& right
		) ;

		std::vector< metro::SampleRange > union_ranges(
			std::vector< metro::SampleRange > const& left,
			metro::SampleRange const& right
		) ;
	}
}

#endif
