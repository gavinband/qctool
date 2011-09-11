#ifndef QCTOOL_CLUSTER_FITTER_HPP
#define QCTOOL_CLUSTER_FITTER_HPP

#include <memory>
#include <vector>
#include <utility>
#include <string>
#include "Eigen/Core"
#include "genfile/SNPDataSourceProcessor.hpp"

struct ClusterFitter: public genfile::SNPDataSourceProcessor::Callback
{
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef std::auto_ptr< ClusterFitter > UniquePtr ;

	static UniquePtr create( appcontext::OptionProcessor const& options, genfile::SNPDataSourceProcessor::Callback* output = 0 ) ;
} ;


#endif
