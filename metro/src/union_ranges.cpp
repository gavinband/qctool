
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <algorithm>
#include "metro/SampleRange.hpp"
#include "metro/union_ranges.hpp"

namespace metro {
	namespace impl {
		std::vector< metro::SampleRange > union_ranges(
			std::vector< metro::SampleRange > const& left,
			std::vector< metro::SampleRange > const& right
		) {
			std::vector< metro::SampleRange > tmp ;
			tmp.reserve( left.size() + right.size() ) ;
			tmp.insert( tmp.end(), left.begin(), left.end() ) ;
			tmp.insert( tmp.end(), right.begin(), right.end() ) ;
			std::sort( tmp.begin(), tmp.end() ) ;
			
			std::vector< metro::SampleRange > result ;
			std::vector< metro::SampleRange >::iterator i = tmp.begin() ;
			while( i != tmp.end() ) {
				std::vector< metro::SampleRange >::iterator rightmost_range = i ;
				std::vector< metro::SampleRange >::iterator j = i ;
				// walk right through the regions that overlap or abut this one.
				for( ++j; j != tmp.end() && j->begin() <= rightmost_range->end(); rightmost_range = j++ ) {
					// nothing to do.
				}
				result.push_back( metro::SampleRange( i->begin(), rightmost_range->end() ) ) ;
				i = j ;
			}
			return result ;
		}

		std::vector< metro::SampleRange > union_ranges(
			std::vector< metro::SampleRange > const& left,
			metro::SampleRange const& right
		) {
			return union_ranges( left, std::vector< metro::SampleRange >( 1, right ) ) ;
		}
	}
}

