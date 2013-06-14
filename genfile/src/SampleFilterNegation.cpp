
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <cassert>
#include "genfile/VariantEntry.hpp"
#include "genfile/SampleFilter.hpp"
#include "genfile/SampleFilterNegation.hpp"

namespace genfile {
	SampleFilterNegation::SampleFilterNegation( SampleFilter::UniquePtr filter ):
		m_filter( filter )
	{
		assert( m_filter.get() ) ;
	}

	void SampleFilterNegation::summarise( std::ostream& o ) const {
		o << "NOT " << ( *m_filter ) ;
	}

	bool SampleFilterNegation::test( genfile::CohortIndividualSource const& source, std::size_t sample, DetailBlock* detail ) const {
		bool result = !m_filter->test( source, sample, detail ) ;
		if( detail ) {
			for( int i = 0; i < detail->cols(); ++i ) {
				(*detail)(0, i) = 1 - (*detail)(0, i) ;
			}
		}
		return result ;
	}
}
