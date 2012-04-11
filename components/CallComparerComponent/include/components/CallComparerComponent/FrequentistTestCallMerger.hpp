
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CALLCOMPARER_COMPONENT_FREQUENTIST_TEST_CALL_MERGER_HPP
#define QCTOOL_CALLCOMPARER_COMPONENT_FREQUENTIST_TEST_CALL_MERGER_HPP

#include <map>
#include <string>
#include <boost/bimap.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <Eigen/Core>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"

struct FrequentistTestCallMerger: PairwiseCallComparerManager::Merger
{
	static SharedPtr create_shared( std::string const& comparison_method, double threshhold ) ;
	FrequentistTestCallMerger( std::string const& comparison_method, double threshhold ) ;
	
	void begin_comparisons( genfile::SNPIdentifyingData const& snp ) ;
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
	std::string m_comparison_method ;
	std::string const m_spec ;
	genfile::SNPIdentifyingData m_snp ;
	typedef std::map< std::pair< std::string, std::string >, genfile::VariantEntry > ComparisonValues ;
	ComparisonValues m_comparison_values ;
	std::pair< double, double > m_range ;
	
	std::set< std::string > m_concordant_calls ;
} ;

#endif

