
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef GENFILE_SAMPLE_FILTER_NEGATION_HPP
#define GENFILE_SAMPLE_FILTER_NEGATION_HPP

#include <memory>
#include <vector>
#include <set>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/VariantEntry.hpp"
#include "SampleFilter.hpp"

namespace genfile {
	struct SampleFilterNegation: public SampleFilter
	{
		SampleFilterNegation( SampleFilter::UniquePtr filter ) ;

		void summarise( std::ostream& ) const ;
		bool test( genfile::CohortIndividualSource const&, std::size_t i, DetailBlock* detail ) const ;
	private:
		SampleFilter::UniquePtr m_filter ;
	} ;
}

#endif

