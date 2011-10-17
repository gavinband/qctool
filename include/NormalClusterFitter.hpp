#ifndef QCTOOL_NORMAL_CLUSTER_FITTER_HPP
#define QCTOOL_NORMAL_CLUSTER_FITTER_HPP

#include <memory>
#include <vector>
#include <utility>
#include <string>
#include "Eigen/Core"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "ClusterFitter.hpp"

struct NormalClusterFitter: public ClusterFitter
{
	NormalClusterFitter( appcontext::OptionProcessor const& options ) ;
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& genotypes ) ;
	void end_processing_snps() ;
private:
	typedef ClusterFitter Base ;
	appcontext::OptionProcessor const& m_options ;
	std::vector< std::pair< std::string, std::string > > const m_spec ;
	double const m_call_threshhold ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
} ;

#endif
