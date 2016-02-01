
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_STRATIFYING_SNP_SUMMARY_COMPUTATION_HPP
#define QCTOOL_STRATIFYING_SNP_SUMMARY_COMPUTATION_HPP

#include <string>
#include <map>
#include <vector>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

struct StratifyingSNPSummaryComputation: public SNPSummaryComputation {
	typedef std::map< genfile::VariantEntry, std::vector< int > > StrataMembers ;
	StratifyingSNPSummaryComputation( SNPSummaryComputation::UniquePtr computation, std::string const& stratification_name, StrataMembers const& strata_members ) ;
	void operator()( VariantIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
private:
	SNPSummaryComputation::UniquePtr m_computation ;
	std::string const m_stratification_name ;
	StrataMembers m_strata_members ;
} ;

#endif
