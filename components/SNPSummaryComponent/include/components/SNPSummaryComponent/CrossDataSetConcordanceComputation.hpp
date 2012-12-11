
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CROSS_DATA_SET_CONCORDANCE_COMPUTATION_HPP
#define QCTOOL_CROSS_DATA_SET_CONCORDANCE_COMPUTATION_HPP

#include <deque>
#include <vector>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSource.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

namespace snp_stats {
	struct CrossDataSetSampleMapper
	{
	public:
		typedef std::vector< genfile::VariantEntry > SampleIdList ;
		typedef std::multimap< std::size_t, std::size_t > SampleMapping ;

		CrossDataSetSampleMapper(
			genfile::CohortIndividualSource const& dataset1_samples,
			std::string const& sample_id_column
		) ;
		
		void set_alternate_dataset(
			genfile::CohortIndividualSource::UniquePtr samples,
			std::string const& comparison_dataset_sample_id_column
		) ;
		
		SampleIdList const& dataset1_sample_ids() const { return m_dataset1_sample_ids ; }
		SampleIdList const& dataset2_sample_ids() const { return m_dataset2_sample_ids ; }
		SampleMapping const& sample_mapping() const { return m_sample_mapping ; }

	private:
		genfile::CohortIndividualSource const& m_dataset1_samples ;
		genfile::CohortIndividualSource::UniquePtr m_dataset2_samples ;
		SampleIdList m_dataset1_sample_ids ;
		SampleIdList m_dataset2_sample_ids ;
		SampleMapping m_sample_mapping ;
	} ;

	struct CrossDataSetComparison: public SNPSummaryComputation
	{
		public:
			typedef std::auto_ptr< CrossDataSetComparison > UniquePtr ;
		
		virtual ~CrossDataSetComparison() {}
		virtual void set_alternate_dataset(
			genfile::CohortIndividualSource::UniquePtr samples,
			std::string const& comparison_dataset_sample_id_column,
			genfile::SNPDataSource::UniquePtr snps
		) = 0 ;
		
		virtual void set_comparer( genfile::SNPIdentifyingData::CompareFields const& comparer ) = 0 ;
	} ;

	struct CrossDataSetConcordanceComputation: public CrossDataSetComparison
	{
	public:
		CrossDataSetConcordanceComputation(
			genfile::CohortIndividualSource const& m_samples,
			std::string const& main_dataset_sample_id_column
		) ;

		void set_alternate_dataset(
			genfile::CohortIndividualSource::UniquePtr samples,
			std::string const& comparison_dataset_sample_id_column,
			genfile::SNPDataSource::UniquePtr snps
		) ;
		
		void set_comparer( genfile::SNPIdentifyingData::CompareFields const& comparer ) ;
		
		void operator()( SNPIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	private:
		CrossDataSetSampleMapper m_sample_mapper ;

		genfile::SNPIdentifyingData::CompareFields m_comparer ;
		genfile::SNPDataSource::UniquePtr m_alt_dataset_snps ;

		double const m_call_threshhold ;
		Eigen::VectorXd m_main_dataset_genotype_subset ;
		Eigen::VectorXd m_alt_dataset_genotypes ;
		Eigen::VectorXd m_pairwise_nonmissingness ;
	} ;

}

#endif
