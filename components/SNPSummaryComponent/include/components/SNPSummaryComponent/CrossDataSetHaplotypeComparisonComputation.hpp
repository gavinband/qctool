
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CROSS_DATA_SET_HAPLOTYPE_COMPARISON_COMPUTATION_HPP
#define QCTOOL_CROSS_DATA_SET_HAPLOTYPE_COMPARISON_COMPUTATION_HPP

#include <deque>
#include <vector>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSource.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/CrossDataSetConcordanceComputation.hpp"

namespace snp_stats {
	struct CrossDataSetHaplotypeComparisonComputation: public CrossDataSetComparison
	{
	public:
		CrossDataSetHaplotypeComparisonComputation(
			genfile::CohortIndividualSource const& m_samples,
			std::string const& main_dataset_sample_id_column
		) ;

		void set_alternate_dataset(
			genfile::CohortIndividualSource::UniquePtr samples,
			std::string const& comparison_dataset_sample_id_column,
			genfile::SNPDataSource::UniquePtr snps
		) ;
	
		void set_comparer( genfile::VariantIdentifyingData::CompareFields const& comparer ) ;
		void set_match_alleles() ;
	
		void operator()( VariantIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	private:
		CrossDataSetSampleMapper m_sample_mapper ;

		genfile::VariantIdentifyingData::CompareFields m_comparer ;
		bool m_match_alleles ;
		genfile::SNPDataSource::UniquePtr m_alt_dataset_snps ;

		double const m_call_threshhold ;
		Eigen::MatrixXd m_haplotypes1 ;
		Eigen::MatrixXd m_haplotypes2 ;
		Eigen::MatrixXd m_nonmissingness1 ;
		Eigen::MatrixXd m_nonmissingness2 ;
	
		Eigen::VectorXd m_relative_phase ;
	} ;
}

#endif
