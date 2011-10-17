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
		virtual void write_intensity_model(
			genfile::SNPIdentifyingData const& snp,
			std::string const& name,
			Eigen::MatrixXd const& intensities,
			genfile::SingleSNPGenotypeProbabilities const& genotypes,
			IntensityModel::SharedPtr model
		) = 0 ;
	} ;

	void connect( ResultCallback::UniquePtr callback ) { m_callbacks.push_back( callback ) ; }
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& ) ;
	void send_results(
		genfile::SNPIdentifyingData const& snp,
		std::string const& name,
		Eigen::MatrixXd const& intensities,
		genfile::SingleSNPGenotypeProbabilities const& genotypes,
		IntensityModel::SharedPtr model
	) ;

private:
	boost::ptr_vector< ResultCallback > m_callbacks ;
} ;

struct IntensityModelSNPCaller {

	
	typedef boost::function< void( std::string const& name, genfile::SingleSNPGenotypeProbabilities const& calls ) > ResultCallback ;

	void connect( ResultCallback )  ;

	void recall_snps(
		genfile::SNPIdentifyingData const& snp,
		std::string const& model_name,
		Eigen::MatrixXd const& intensities,
		genfile::SingleSNPGenotypeProbabilities const& genotypes,
		IntensityModel::SharedPtr model
	) ;

private:
	boost::ptr_vector< ResultCallback > m_callbacks ;
} ;

#endif
