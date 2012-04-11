
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_COMPOSITE_COHORT_INDIVIDUAL_SOURCE_HPP
#define GENFILE_COMPOSITE_COHORT_INDIVIDUAL_SOURCE_HPP

#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	// This class is a base class for cohort individual sources which form a hierarchical view of the represented data.
	class CompositeCohortIndividualSource: public CohortIndividualSource
	{
		// Find the parent of this source
		virtual CohortIndividualSource const& get_parent_source() const = 0;
	} ;
}

#endif
