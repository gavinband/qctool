#ifndef QCTOOL_CLUSTER_FITTER_HPP
#define QCTOOL_CLUSTER_FITTER_HPP

#include <memory>
#include <vector>
#include <utility>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include "Eigen/Core"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "IntensityModel.hpp"

struct ClusterFitter: public genfile::SNPDataSourceProcessor::Callback
{
public:
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef std::auto_ptr< ClusterFitter > UniquePtr ;

	static UniquePtr create( appcontext::OptionProcessor const& options ) ;

	struct ResultCallback {
		typedef std::auto_ptr< ResultCallback > UniquePtr ;
		virtual ~ResultCallback() {}
		virtual void set_SNP( genfile::SNPIdentifyingData const& ) = 0 ;
		virtual void write_cluster_fit( std::string const& name, std::vector< std::size_t > const&, Eigen::MatrixXd const& fit ) = 0 ;
	} ;

	void connect( ResultCallback::UniquePtr callback ) { m_callbacks.push_back( callback ) ; }
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& ) ;
	void send_results( std::string const& name, std::vector< std::size_t > const& non_missing_counts, Eigen::MatrixXd const& fit ) ;

private:
	boost::ptr_vector< ResultCallback > m_callbacks ;
} ;


#endif
