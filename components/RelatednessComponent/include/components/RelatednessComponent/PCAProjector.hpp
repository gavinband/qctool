#ifndef RELATEDNESS_COMPONENT_PCA_PROJECTOR_HPP
#define RELATEDNESS_COMPONENT_PCA_PROJECTOR_HPP

#include <map>
#include <boost/signals2/signal.hpp>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "appcontext/UIContext.hpp"

namespace pca {
	struct PCAProjector: public genfile::SNPDataSourceProcessor::Callback
	{
	public:
		typedef Eigen::MatrixXd Matrix ;
		typedef Eigen::VectorXd Vector ;
		typedef std::auto_ptr< PCAProjector > UniquePtr ;
	public:
		PCAProjector( genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) ;
		
		void set_loadings( std::vector< genfile::SNPIdentifyingData > const& snps, Matrix const& loadings, std::vector< std::string > const& names ) ;

		void begin_processing_snps( std::size_t number_of_samples ) ;
		void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& ) ;
		void end_processing_snps() ;

		typedef boost::function< genfile::VariantEntry ( std::size_t ) > GetNames ;
		typedef boost::signals2::signal< void( std::string, Eigen::VectorXd const&, Eigen::MatrixXd const&, GetNames, GetNames ) > ResultSignal ;
		void send_results_to( ResultSignal::slot_type callback ) ;
		void send_results( std::string description, Eigen::MatrixXd const& projections, GetNames pca_row_names, GetNames pca_column_names ) ;

	private:
		genfile::CohortIndividualSource const& m_samples ;
		appcontext::UIContext& m_ui_context ;
		typedef std::map< genfile::SNPIdentifyingData, int > SnpMap ;
		SnpMap m_snps ;
		Eigen::MatrixXd m_loadings ;
		std::vector< std::string > m_names ;
		Eigen::VectorXd m_genotype_calls ;
		Eigen::VectorXd m_non_missingness ;
		Eigen::MatrixXd m_projections ;
	
		// m_visited keeps track of which SNPs have been used in the projection.
		typedef std::map< genfile::SNPIdentifyingData, bool > VisitedSnpMap ;
		VisitedSnpMap m_visited ;
	
		ResultSignal m_result_signal ;
		
		void diagnose_projection() const ;
	} ;
}

#endif

