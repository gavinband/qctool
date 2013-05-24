
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef GENFILE_VARIABLEINRANGE_SAMPLE_FILTER_HPP
#define GENFILE_VARIABLEINRANGE_SAMPLE_FILTER_HPP

#include <memory>
#include <vector>
#include <set>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SampleFilter.hpp"
#include "genfile/impl/cast_types_for_comparison.hpp"

namespace genfile {
	struct VariableInRangeSampleFilter: public SampleFilter
	{
		VariableInRangeSampleFilter( std::string const& statistic_name, genfile::VariantEntry lower_bound, genfile::VariantEntry upper_bound ) ;

		genfile::VariantEntry const& lower_bound() const { return m_lower_bound ; }
		genfile::VariantEntry const& upper_bound() const { return m_upper_bound ; }
		std::string const& variable() const { return m_variable ; }

	private:
		std::string const m_variable ;
		genfile::VariantEntry m_lower_bound ;
		genfile::VariantEntry m_upper_bound ;
	
		bool test( genfile::CohortIndividualSource const& source, std::size_t i ) const ;
		virtual bool compare_impl( genfile::VariantEntry const& value, genfile::VariantEntry const& lower, genfile::VariantEntry const& upper ) const = 0 ;
	} ;

	struct VariableInInclusiveRangeSampleFilter: public VariableInRangeSampleFilter
	{
		VariableInInclusiveRangeSampleFilter( std::string const& variable, genfile::VariantEntry lower_bound, genfile::VariantEntry upper_bound ) ;
		void summarise( std::ostream& oStream ) const ;
	private:
		bool compare_impl( genfile::VariantEntry const& value, genfile::VariantEntry const& lower, genfile::VariantEntry const& upper ) const ;
	} ;

	struct VariableInExclusiveRangeSampleFilter: public VariableInRangeSampleFilter
	{
		VariableInExclusiveRangeSampleFilter( std::string const& variable, genfile::VariantEntry lower_bound, genfile::VariantEntry upper_bound ) ;
		void summarise( std::ostream& oStream ) const ;
	private:
		bool compare_impl( genfile::VariantEntry const& value, genfile::VariantEntry const& lower, genfile::VariantEntry const& upper ) const ;
	} ;
}

#endif

