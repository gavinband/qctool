
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef GENFILE_IMPL_CAST_TYPES_FOR_COMPARISON_HPP
#define GENFILE_IMPL_CAST_TYPES_FOR_COMPARISON_HPP

#include "genfile/VariantEntry.hpp"
#include "genfile/string_utils/string_utils.hpp"

namespace genfile {
	namespace impl {
		// 
		// Homogenise types for value comparison.
		// We do this:
		// * 1. convert any ints to doubles
		// * 2. of one operand is int or double and the other operand is string then we cast the string to a double if possible.
		//
		void cast_types_for_comparison( genfile::VariantEntry& value1, genfile::VariantEntry& value2 ) ;
	}
}

#endif

