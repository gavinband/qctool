
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGEN_REFERENCE_IMPLEMENTATION_TYPES_HPP
#define BGEN_REFERENCE_IMPLEMENTATION_TYPES_HPP

#include <stdint.h>
#include <exception>
#include "../../config.hpp"

/*
* This file contains a reference implementation of the BGEN file format
* specification described at:
* http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format.html
*
* To use this file you will also need "endianness_utils.hpp".
*
*/

namespace genfile {
	namespace bgen {
		struct BGenError: public virtual std::exception {
			~BGenError() throw() {}
			char const* what() const throw() { return "BGenError" ; }
		} ;
		typedef ::uint32_t uint32_t ;
		typedef ::uint16_t uint16_t ;

		// v1.1 definitions were:
		// enum FlagsType { e_NoFlags = 0, e_CompressedSNPBlocks = 0x1, e_LongIds = 0x4 } ;
		// v1.2 definitions:
		enum FlagMask { e_NoFlags = 0, e_CompressedSNPBlocks = 0x1, e_Layout = 0x3C } ;
		enum Layout { e_v10Layout = 0x0, e_v11Layout = 0x4, e_v12Layout = 0x8 } ;
		enum Structure { e_SampleIdentifiers = 0x80000000 } ;
	}
}

#endif
