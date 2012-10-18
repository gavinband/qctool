
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
	struct CrossDataSetConcordanceComputation: public SNPSummaryComputation
	{
	public:
		typedef std::auto_ptr< CrossDataSetConcordanceComputation > UniquePtr ;
		typedef std::vector< genfile::VariantEntry > SampleIdList ;
		typedef std::multimap< std::size_t, std::size_t > SampleMapping ;
		
	public:
		CrossDataSetConcordanceComputation(
			genfile::CohortIndividualSource const& m_samples
		) ;

		void set_alternate_dataset( genfile::CohortIndividualSource::UniquePtr samples, genfile::SNPDataSource::UniquePtr snps ) ;
		void operator()( SNPIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	private:
		genfile::CohortIndividualSource const& m_samples ;
		SampleIdList m_sample_ids ;
		double const m_call_threshhold ;
		genfile::CohortIndividualSource::UniquePtr m_alt_dataset_samples ;
		genfile::SNPDataSource::UniquePtr m_alt_dataset_snps ;
		SampleIdList m_alt_dataset_sample_ids ;
		SampleMapping m_sample_mapping ;

		Eigen::VectorXd m_alt_dataset_genotypes ;
		Eigen::VectorXd m_concordance ;
	} ;
}

#endif
