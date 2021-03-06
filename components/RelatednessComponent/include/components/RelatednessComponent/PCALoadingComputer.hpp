
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
	void set_UDUT( std::size_t number_of_snps, Matrix const& udut_decomposition ) ;
	void set_number_of_loadings( std::size_t n ) ;

	void begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) ;
	void processed_snp( genfile::VariantIdentifyingData const&, genfile::VariantDataReader& ) ;
	void end_processing_snps() ;

	typedef boost::function< genfile::VariantEntry ( std::size_t ) > GetNames ;
	typedef boost::signals2::signal< void( genfile::VariantIdentifyingData const&, double const, double const, Eigen::VectorXd const&, GetNames ) > ResultSignal ;
	typedef ResultSignal::slot_type ResultCallback ;

	void send_results_to( ResultCallback callback ) ;
	void send_results( genfile::VariantIdentifyingData const& snp, double const, double const, Eigen::VectorXd const& data, GetNames ) ;
	
	std::string get_metadata() const ;
private:
	Eigen::VectorXd m_D ;
	Eigen::VectorXd m_sqrt_D_inverse ;
	Eigen::MatrixXd m_U ;
	int const m_number_of_loadings ;
	int m_number_of_snps ;
	Eigen::RowVectorXd m_loading_vectors ;
	Eigen::VectorXd m_allele_frequencies ;
	Eigen::VectorXd m_genotype_calls ;
	Eigen::VectorXd m_non_missingness ;
	
	ResultSignal m_result_signal ;
} ;

#endif

