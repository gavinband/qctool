
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_RISK_SCORE_COMPUTATION_HPP
#define QCTOOL_RISK_SCORE_COMPUTATION_HPP

#include <string>
#include <map>
#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <Eigen/Core>
#include "metro/mean_and_variance.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace sample_stats {
	struct RiskScoreComputation: public SampleSummaryComputation
	{
		typedef std::auto_ptr< RiskScoreComputation > UniquePtr ;
		RiskScoreComputation( genfile::CohortIndividualSource const& samples, genfile::SNPIdentifyingData::CompareFields comparator ) ;
		void accumulate( genfile::SNPIdentifyingData const&, Genotypes const&, genfile::VariantDataReader& ) ;
		void compute( ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

		typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;
		void add_effects( statfile::BuiltInTypeStatSource& source, ProgressCallback ) ;
		
	private:
		std::size_t const m_number_of_samples ;

		typedef std::set< std::string > Identifiers ;
		typedef Eigen::MatrixXd Betas ;
		typedef std::map< std::string, Betas > RiskScoreIdBetaMap ;
		typedef std::map< genfile::SNPIdentifyingData, RiskScoreIdBetaMap, genfile::SNPIdentifyingData::CompareFields > RiskScoreMap ;
		
		Identifiers m_identifiers ;
		RiskScoreMap m_map ;

		std::map< std::string, std::size_t > m_snp_counts ;
		std::map< std::string, Eigen::VectorXd > m_scores ;
		std::map< std::string, Eigen::VectorXd > m_counts ;
	private:
		
		void accumulate( Genotypes const& genotypes, Betas const& betas ) ;
		void accumulate( Genotypes const& genotypes, Betas const& betas, Eigen::VectorXd& scores, Eigen::VectorXd& counts ) const ;
		void compute_impl( std::string const& risk_score_identifier, ResultCallback callback ) ;
	} ;
}

#endif
