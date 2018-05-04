
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_SAMPLEDATAINTERPRETATION_HPP
#define METRO_SAMPLEDATAINTERPRETATION_HPP

#include <string>
#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>
#include <Eigen/Core>
#include "metro/SampleRange.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CrossCohortCovariateValueMapping.hpp"

namespace metro {
	/* The purpose of this struct is to
	* convert data in sample file representation
	* into a matrix of values suitable for regression
	* This means laying out the data in a matrix
	* and expanding discrete columns to multiple columns.
	*/
	struct SampleDataInterpretation {

		typedef Eigen::MatrixXd Matrix ;
		typedef std::vector< metro::SampleRange > SampleRanges ;

		SampleDataInterpretation( genfile::CohortIndividualSource const& samples, std::vector< std::string > column_names ) ;
		void add_samples( genfile::CohortIndividualSource const& samples ) ;

		Matrix compute_result() const ;
		std::vector< std::string > const& column_names() const ;
		SampleRanges compute_nonmissing_samples() const ;
	private:
		genfile::CohortIndividualSource const& m_samples ;
		genfile::CohortIndividualSource::ColumnSpec const m_column_spec ;
		std::vector< std::string > const m_columns ;
		boost::ptr_vector< genfile::CrossCohortCovariateValueMapping > m_mappings ;
	} ;
}

#endif
