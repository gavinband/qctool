
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CALLCOMPARER_COMPONENT_ACCEPT_ALL_CALL_MERGER_HPP
#define QCTOOL_CALLCOMPARER_COMPONENT_ACCEPT_ALL_CALL_MERGER_HPP

#include <map>
#include <string>
#include <boost/bimap.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <Eigen/Core>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "components/SNPSummaryComponent/PairwiseCallComparerManager.hpp"

// A call merger that accepts all calls passed to it.
struct AcceptAllCallMerger: PairwiseCallComparerManager::Merger
{
	AcceptAllCallMerger() ;
	
	void begin_comparisons( genfile::SNPIdentifyingData const& snp ) ;
	void add_callset( std::string const& ) ;
	void set_result(
		std::string const& callset1,
		std::string const& callset2,
		std::string const& comparison_method,
		std::string const& comparison_variable,
		genfile::VariantEntry const& value
	) ;
	void end_comparisons() ;

	std::string get_spec() const ;
	std::string get_result_as_string() const ;
private:
	std::set< std::string > m_calls ;
} ;

#endif

