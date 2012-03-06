#ifndef RELATEDNESS_COMPONENT_PCA_LOADING_COMPUTER_HPP
#define RELATEDNESS_COMPONENT_PCA_LOADING_COMPUTER_HPP

#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "components/RelatednessComponent/KinshipCoefficientManager.hpp"

struct PCALoadingComputer: public genfile::SNPDataSourceProcessor::Callback
{
public:
	typedef Eigen::MatrixXd Matrix ;
	typedef std::auto_ptr< PCALoadingComputer > UniquePtr ;
public:
	// set the PCA components.  This should be a Nxl matrix
	// where N is the number of samples and l the number of PCA components.
	void set_PCA_components( Matrix const& pca_components ) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& ) ;
	void end_processing_snps() ;

	typedef boost::function< genfile::VariantEntry ( std::size_t ) > GetNames ;
	typedef boost::function< void( genfile::SNPIdentifyingData const&, Eigen::VectorXd const&, GetNames ) > ResultCallback ;
	void send_results_to( ResultCallback callback ) ;
	void send_results( genfile::SNPIdentifyingData const& snp, Eigen::VectorXd const& data, GetNames ) ;
	
private:	
	Eigen::MatrixXd m_PCA_components ;
	Eigen::VectorXd m_loading_vectors ;
	Eigen::VectorXd m_genotype_calls ;
	Eigen::VectorXd m_non_missingness ;
	
	typedef boost::signals2::signal< void( genfile::SNPIdentifyingData const&, Eigen::VectorXd const&, GetNames ) > ResultSignal ;
	ResultSignal m_result_signal ;
} ;

#endif

