
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <map>
#include <string>
#include <boost/bimap.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <Eigen/Core>

#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "components/SNPSummaryComponent/PairwiseCallComparerManager.hpp"
#include "components/SNPSummaryComponent/AcceptAllCallMerger.hpp"

AcceptAllCallMerger::AcceptAllCallMerger()
{}

std::string AcceptAllCallMerger::get_spec() const {
	return "AcceptAllCallMerger" ;
}

std::string AcceptAllCallMerger::get_result_as_string() const {
	std::string result ;
	foreach( std::string const& call, m_calls ) {
		if( result.size() > 0 ) {
			result += "," ;
		}
		result += call ;
	}
	return result ;
}

void AcceptAllCallMerger::begin_comparisons( genfile::VariantIdentifyingData const& ) {
	m_calls.clear() ;
}

void AcceptAllCallMerger::add_callset( std::string const&  callset ) {
	m_calls.insert( callset ) ;
}

void AcceptAllCallMerger::set_result(
	std::string const& callset1,
	std::string const& callset2,
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
)
{}

void AcceptAllCallMerger::end_comparisons() {
}
