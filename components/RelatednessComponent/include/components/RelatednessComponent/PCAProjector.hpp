
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

		static genfile::VariantEntry get_projection_name( std::size_t i ) ;

	public:
		static UniquePtr create(
			genfile::CohortIndividualSource const& samples,
			appcontext::UIContext& ui_context,
			genfile::VariantIdentifyingData::CompareFields const&
		) ;
		PCAProjector(
			genfile::CohortIndividualSource const& samples,
			appcontext::UIContext& ui_context,
			genfile::VariantIdentifyingData::CompareFields const&
		) ;
		
		void set_loadings(
			std::vector< genfile::VariantIdentifyingData > const& snps,
			Vector const& counts,
			Vector const& frequencies,
			Matrix const& loadings,
			std::vector< std::string > const& names
		) ;

		void begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) ;
		void processed_snp( genfile::VariantIdentifyingData const&, genfile::VariantDataReader& ) ;
		void end_processing_snps() ;

		typedef boost::function< genfile::VariantEntry ( std::size_t ) > GetNames ;
		typedef boost::signals2::signal< void( std::string, Eigen::MatrixXd const&, GetNames, GetNames ) > ResultSignal ;
		void send_results_to( ResultSignal::slot_type callback ) ;

		std::string get_metadata() const ;
	private:
		genfile::CohortIndividualSource const& m_samples ;
		appcontext::UIContext& m_ui_context ;
		typedef std::map< genfile::VariantIdentifyingData, int, genfile::VariantIdentifyingData::CompareFields > SnpMap ;
		SnpMap m_snps ;
		Eigen::VectorXd m_counts ;
		Eigen::VectorXd m_frequencies ;
		Eigen::MatrixXd m_loadings ;
		std::vector< std::string > m_names ;
		Eigen::VectorXd m_genotype_calls ;
		Eigen::VectorXd m_non_missingness ;
		Eigen::MatrixXd m_projections ;
		Eigen::VectorXd m_snps_visited_per_sample ;
	
		// m_visited keeps track of which SNPs have been used in the projection.
		typedef std::map< genfile::GenomePosition, int > VisitedSnpMap ;
		VisitedSnpMap m_visited ;
		std::size_t m_total_snps_visited ;
		ResultSignal m_result_signal ;
		
		void diagnose_projection() const ;
	} ;
}

#endif

