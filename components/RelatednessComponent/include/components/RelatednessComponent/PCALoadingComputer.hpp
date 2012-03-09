#ifndef RELATEDNESS_COMPONENT_PCA_LOADING_COMPUTER_HPP
#define RELATEDNESS_COMPONENT_PCA_LOADING_COMPUTER_HPP

#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "components/RelatednessComponent/KinshipCoefficientManager.hpp"

struct PCALoadingComputer: public genfile::SNPDataSourceProcessor::Callback
{
public:
	typedef Eigen::MatrixXd Matrix ;
	typedef Eigen::VectorXd Vector ;
	typedef std::auto_ptr< PCALoadingComputer > UniquePtr ;
public:
	PCALoadingComputer( int number_of_loadings ) ;
	void set_UDUT( Matrix const& udut_decomposition ) ;
	void set_number_of_loadings( std::size_t n ) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& ) ;
	void end_processing_snps() ;

	typedef boost::function< genfile::VariantEntry ( std::size_t ) > GetNames ;
	typedef boost::function< void( genfile::SNPIdentifyingData const&, Eigen::VectorXd const&, GetNames ) > ResultCallback ;
	void send_results_to( ResultCallback callback ) ;
	void send_results( genfile::SNPIdentifyingData const& snp, Eigen::VectorXd const& data, GetNames ) ;
	
private:
	Eigen::MatrixXd m_D ;
	Eigen::MatrixXd m_U ;
	int const m_number_of_loadings ;
	int m_number_of_snps ;
	Eigen::VectorXd m_loading_vectors ;
	Eigen::VectorXd m_genotype_calls ;
	Eigen::VectorXd m_non_missingness ;
	
	typedef boost::signals2::signal< void( genfile::SNPIdentifyingData const&, Eigen::VectorXd const&, GetNames ) > ResultSignal ;
	ResultSignal m_result_signal ;
} ;

#endif

