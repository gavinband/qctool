
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_INTENSITY_REPORTER_HPP
#define QCTOOL_INTENSITY_REPORTER_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantDataReader.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

namespace snp_summary_component {
	struct IntensityReporter: public SNPSummaryComputation {
		IntensityReporter( genfile::CohortIndividualSource const& samples, double call_threshhold = 0.9 ) ;
		void operator()( SNPIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
	private:
		genfile::CohortIndividualSource const& m_samples ;
		std::vector< genfile::VariantEntry > m_sample_ids ;
		double const m_call_threshhold ;
		typedef Eigen::MatrixXd IntensityMatrix ;
		IntensityMatrix m_intensities ;
		IntensityMatrix m_nonmissingness ;
	} ;
}

#endif
