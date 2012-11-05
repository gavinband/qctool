
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_LEASTMISSINGCALLCOMPARERCOMPONENT_CONSENSUS_CALLER_HPP
#define QCTOOL_LEASTMISSINGCALLCOMPARERCOMPONENT_CONSENSUS_CALLER_HPP

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "components/SNPSummaryComponent/ConsensusCaller.hpp"

// Given a set of calls passed in by processed_snp
// The consensus caller outputs a 
struct LeastMissingConsensusCaller: public ConsensusCaller
{
public:
	LeastMissingConsensusCaller() ;
public:
	void set_result(
		std::string const& comparison,
		std::string const& accepted_calls,
		PairwiseCallComparerManager::Calls const& calls
	) ;

private:
	double const m_call_threshhold ;
	std::vector< Eigen::MatrixXd > m_genotypes ;
} ;

#endif
