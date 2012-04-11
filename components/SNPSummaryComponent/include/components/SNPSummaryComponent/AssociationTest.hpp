
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_ASSOCIATION_TEST_HPP
#define QCTOOL_ASSOCIATION_TEST_HPP

#include <vector>
#include <string>
#include <limits>
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

struct AssociationTest: public SNPSummaryComputation {
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;

	static UniquePtr create(
		std::string const& phenotype_name,
		std::vector< std::string > const& covariate_names,
		genfile::CohortIndividualSource const& samples,
		appcontext::OptionProcessor const& options
	) ;

private:
} ;

#endif
