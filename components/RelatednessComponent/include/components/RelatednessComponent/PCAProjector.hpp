
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef RELATEDNESS_COMPONENT_PCA_PROJECTOR_HPP
#define RELATEDNESS_COMPONENT_PCA_PROJECTOR_HPP

#include <map>
#include <memory>
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
		static UniquePtr create( genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) ;
		PCAProjector( genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) ;
		
		void set_loadings( std::vector< genfile::SNPIdentifyingData > const& snps, Matrix const& loadings, std::vector< std::string > const& names ) ;

		void begin_processing_snps( std::size_t number_of_samples ) ;
		void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& ) ;
		void end_processing_snps() ;

		typedef boost::function< genfile::VariantEntry ( std::size_t ) > GetNames ;
		typedef boost::signals2::signal< void( std::string, Eigen::MatrixXd const&, GetNames, GetNames ) > ResultSignal ;
		void send_results_to( ResultSignal::slot_type callback ) ;

		std::string get_metadata() const ;
	private:
		genfile::CohortIndividualSource const& m_samples ;
		appcontext::UIContext& m_ui_context ;
		typedef std::map< genfile::GenomePosition, int > SnpMap ;
		SnpMap m_snps ;
		Eigen::MatrixXd m_loadings ;
		std::vector< std::string > m_names ;
		Eigen::VectorXd m_genotype_calls ;
		Eigen::VectorXd m_non_missingness ;
		Eigen::MatrixXd m_projections ;
	
		// m_visited keeps track of which SNPs have been used in the projection.
		typedef std::map< genfile::GenomePosition, bool > VisitedSnpMap ;
		VisitedSnpMap m_visited ;
	
		ResultSignal m_result_signal ;
		
		void diagnose_projection() const ;
	} ;
}

#endif
