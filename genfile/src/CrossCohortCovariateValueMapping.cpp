
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CrossCohortCovariateValueMapping.hpp"
#include "genfile/NormalisingCrossCohortCovariateValueMapping.hpp"
#include "genfile/CategoricalCrossCohortCovariateValueMapping.hpp"
#include "genfile/ContinuousVariableCrossCohortCovariateValueMapping.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	CrossCohortCovariateValueMapping::UniquePtr CrossCohortCovariateValueMapping::create( CohortIndividualSource::SingleColumnSpec const& column_spec ) {
		if( column_spec.is_continuous() ) {
			if( column_spec.is_phenotype() ) {
				return UniquePtr( new NormalisingCrossCohortCovariateValueMapping( column_spec.name() )) ;
			}
			else {
				return UniquePtr( new ContinuousVariableCrossCohortCovariateValueMapping( column_spec.name() )) ;
			}
		}
		else {
			return UniquePtr( new CategoricalCrossCohortCovariateValueMapping( column_spec.name() )) ;
		}
	}
}
