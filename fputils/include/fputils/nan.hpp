
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef FPUTILS_NAN_HPP
#define FPUTILS_NAN_HPP

namespace fputils {
	bool is_NaN( double value ) {
		return value != value ;
	}
}

#endif
