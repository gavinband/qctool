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
