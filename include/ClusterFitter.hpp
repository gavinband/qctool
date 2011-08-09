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

struct NormalClusterFitter: public ClusterFitter
{
	NormalClusterFitter( appcontext::OptionProcessor const& options ) ;
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& genotypes ) ;
	void end_processing_snps() ;
private:
	appcontext::OptionProcessor const& m_options ;
	std::string const m_filename ;
	std::auto_ptr< std::ostream > const m_file ;
	std::vector< std::pair< std::string, std::string > > const m_spec ;
	double const m_call_threshhold ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;

	void get_cluster_fit(
		std::vector< genfile::VariantEntry > const& genotypes,
		Eigen::MatrixXd const& intensities,
		std::vector< Eigen::Vector2d >& means,
		std::vector< Eigen::Matrix2d >& variances,
		std::vector< std::size_t >& non_missing_counts
	) ;

	void write_output_columns(
		std::vector< std::pair< std::string, std::string > > const& spec
	) const ;

	void write_id_data( genfile::SNPIdentifyingData const& id_data ) const ;

	void write_cluster_fit(
		std::vector< Eigen::Vector2d > const& means,
		std::vector< Eigen::Matrix2d > const& variances,
		std::vector< std::size_t > const& non_missing_counts
	) const ;
} ;

#endif
