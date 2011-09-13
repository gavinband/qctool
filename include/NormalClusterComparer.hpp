#ifndef QCTOOL_NORMAL_CLUSTER_FIT_COMPARER_HPP
#define QCTOOL_NORMAL_CLUSTER_FIT_COMPARER_HPP

#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "ClusterFitter.hpp"

struct NormalClusterFitComparer: public boost::noncopyable
{
	typedef std::auto_ptr< NormalClusterFitComparer > UniquePtr ;
	virtual ~NormalClusterFitComparer() {}
	virtual void compare( Eigen::MatrixXd const& intensities, Eigen::MatrixXd const& fit1, Eigen::MatrixXd const& fit2 ) = 0 ;
} ;


struct NormalClusterFitComparerThing: public ClusterFitter:ResultCallback, public genfile::SNPDataSourceProcessor::Callback {
public:
	virtual ~NormalClusterFitComparerThing() {} ;

	void add_comparer( std::string const& name, NormalClusterFitComparer::UniquePtr comparer ) ;
	void processed_snp( genfile::SNPIdentifyingData snp, genfile::VariantDataReader& data_reader ) ;
	
private:
	boost::ptr_map< std::string, NormalClusterFitComparer > m_comparers ;
} ;

#endif
