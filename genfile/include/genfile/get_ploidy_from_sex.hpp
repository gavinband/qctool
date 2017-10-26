
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_GET_PLOIDY_FROM_SEX_HPP
#define GENFILE_GET_PLOIDY_FROM_SEX_HPP

#include <vector>
#include "genfile/VariantIdentifyingData.hpp"

namespace genfile {
	// Compute ploidy from the sex of a sample and chromosome.
	int get_ploidy_from_sex(
		std::vector< char > const& sexes,
		genfile::VariantIdentifyingData const& variant,
		std::size_t sample_i
	) ;
}

#endif
