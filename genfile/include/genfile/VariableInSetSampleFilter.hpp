
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef GENFILE_VARIABLEINSET_SAMPLE_FILTER_HPP
#define GENFILE_VARIABLEINSET_SAMPLE_FILTER_HPP

#include <memory>
#include <vector>
#include <set>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SampleFilter.hpp"
#include "genfile/impl/cast_types_for_comparison.hpp"

namespace genfile {
	struct VariableInSetSampleFilter: public SampleFilter
	{
	public:
		typedef std::auto_ptr< VariableInSetSampleFilter > UniquePtr ;
		
	public:
		VariableInSetSampleFilter( std::string const& variable ) ;

		void summarise( std::ostream& o ) const ;
		void add_level( genfile::VariantEntry value ) ;
		bool test( genfile::CohortIndividualSource const& source, std::size_t i, DetailBlock* ) const ;

	private:
		std::string const m_variable ;
		std::vector< genfile::VariantEntry > m_levels ;
	} ;
}

#endif
