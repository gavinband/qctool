
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include "metro/UnivariateLogDensity.hpp"
#include "metro/distributions/Flat.hpp"

namespace metro {
	namespace distributions {
		Flat::UniquePtr create() {
			return Flat::UniquePtr( new Flat ) ;
		}
		
		std::string Flat::get_summary() const {
			return "flat" ;
		}
	}
}

