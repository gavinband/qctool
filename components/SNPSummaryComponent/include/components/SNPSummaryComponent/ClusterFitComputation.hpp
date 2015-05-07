
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CLUSTER_FIT_COMPUTATION_HPP
#define QCTOOL_CLUSTER_FIT_COMPUTATION_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantDataReader.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

namespace snp_summary_component {
	struct ClusterFitComputation: public SNPSummaryComputation {
		typedef std::auto_ptr< ClusterFitComputation > UniquePtr ;
		ClusterFitComputation(
			double nu,
			double regularisationVariance = 0.5,
			double regularisationWeight = 10,
			double call_threshhold = 0.9
		) ;
		void operator()( SNPIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
		void set_scale( std::string const& scale ) ;
	private:
		double const m_call_threshhold ;
		double const m_nu ;
		std::string m_scale ;
		std::string m_xAxisName ;
		std::string m_yAxisName ;
		
		Eigen::MatrixXd const m_regularisingSigma ;
		double const m_regularisingWeight ;
		typedef Eigen::MatrixXd IntensityMatrix ;
		IntensityMatrix m_intensities ;
		IntensityMatrix m_nonmissingness ;
	} ;
}

#endif
