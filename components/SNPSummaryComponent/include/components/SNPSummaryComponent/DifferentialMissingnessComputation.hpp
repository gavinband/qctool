
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_DIFFERENTIAL_MISSINGNESS_COMPUTATION_HPP
#define QCTOOL_DIFFERENTIAL_MISSINGNESS_COMPUTATION_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantDataReader.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

struct DifferentialMissingnessComputation: public SNPSummaryComputation {
	typedef std::map< genfile::VariantEntry, std::vector< int > > StrataMembers ;
	static UniquePtr create( std::string const& stratification_name, StrataMembers const& strata_members ) ;
	DifferentialMissingnessComputation( std::string const& stratification_name, StrataMembers const& strata_members, double threshhold = 0.9 ) ;
	void operator()( SNPIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
private:
	std::string const m_stratification_name ;
	StrataMembers const m_strata_members ;
	std::vector< int > const m_strata_levels ;
	double const m_threshhold ;
private:
	std::vector< int > compute_strata_levels( StrataMembers const& strata_members ) const ;
} ;

#endif
