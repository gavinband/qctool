
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_QUANGSTYLECALLCOMPARERCOMPONENT_CONSENSUS_CALLER_HPP
#define QCTOOL_QUANGSTYLECALLCOMPARERCOMPONENT_CONSENSUS_CALLER_HPP

#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "components/SNPSummaryComponent/ConsensusCaller.hpp"

// Given a set of calls passed in by processed_snp
// The consensus caller outputs a 
struct QuangStyleConsensusCaller: public ConsensusCaller
{
public:
	QuangStyleConsensusCaller( double threshhold, unsigned int minimum_consensus ) ;
public:
	void set_result(
		std::string const& comparison,
		std::string const& accepted_calls,
		PairwiseCallComparerManager::Calls const& calls
	) ; 
private:
	double const m_call_threshhold ;
	double const m_minimum_consensus ;
	Eigen::VectorXd m_genotypes ;
	Eigen::VectorXd m_result_calls ;
	Eigen::VectorXd m_consensus_counts ;
	Eigen::MatrixXd m_result_probs ;
} ;

#endif
