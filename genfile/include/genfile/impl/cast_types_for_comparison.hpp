
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
		void cast_types_for_comparison( genfile::VariantEntry& value1, genfile::VariantEntry& value2 ) {
			using genfile::string_utils::to_repr ;
			using genfile::string_utils::StringConversionError ;
			if( value1.is_int() ) {
				value1 = double( value1.as< int >() ) ;
			}
			if( value2.is_int() ) {
				value2 = double( value2.as< int >() ) ;
			}

			if( value1.is_double() && value2.is_string() ) {
				try {
					value2 = to_repr< double >( value2.as< std::string >() ) ;
				} catch( StringConversionError const& ) {
				}
			}
			else if( value2.is_double() && value1.is_string() ) {
				try {
					value1 = to_repr< double >( value1.as< std::string >() ) ;
				} catch( StringConversionError const& ) {
				}
			}
		}
	}
}

#endif

