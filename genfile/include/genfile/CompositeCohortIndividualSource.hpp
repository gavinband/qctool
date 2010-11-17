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
