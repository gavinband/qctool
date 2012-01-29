#ifndef QCTOOL_ASSOCIATION_TEST_HPP
#define QCTOOL_ASSOCIATION_TEST_HPP

#include <vector>
#include <string>
#include <limits>
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "SNPSummaryComputation.hpp"

struct AssociationTest: public SNPSummaryComputation {
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;

	typedef std::auto_ptr< AssociationTest > UniquePtr ;
	static UniquePtr create(
		std::string const& phenotype_name,
		std::vector< std::string > const& covariate_names,
		genfile::CohortIndividualSource const& samples,
		appcontext::OptionProcessor const& options
	) ;

private:
	static Vector get_sample_column(
		genfile::CohortIndividualSource const& samples,
		std::string const& column_name
	) {
		Vector result = Vector::Constant( samples.get_number_of_individuals(), std::numeric_limits< double >::quiet_NaN() ) ;
		for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
			genfile::CohortIndividualSource::Entry const& entry = samples.get_entry( i, column_name ) ;
			if( !entry.is_missing() ) {
				result(i) = entry.as< double >() ;
			}
		}
		return result ;
	}
} ;

#endif
