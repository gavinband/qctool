
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "metro/regression/LogLikelihood.hpp"

namespace metro {
	namespace regression {
		LogLikelihood::~LogLikelihood() {}
		
		std::vector< std::string > LogLikelihood::get_parameter_names() const {
			std::size_t const n = identify_parameters().rows() ;
			std::vector< std::string > result(n, "") ;
			for( std::size_t i = 0; i < n; ++i ) {
				result[i] = get_parameter_name(i) ;
			}
			return result ;
		}
	}
}

