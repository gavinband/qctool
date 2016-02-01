
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPPOSITIONINRANGETEST_HPP
#define SNPPOSITIONINRANGETEST_HPP

#include <set>
#include <string>
#include "genfile/GenomePosition.hpp"
#include "genfile/GenomePositionRange.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"

namespace genfile {
	// Test if the given SNP is in the given (closed) range of positions.
	struct SNPPositionInRangeTest: public VariantIdentifyingDataTest
	{
		SNPPositionInRangeTest( GenomePositionRange const range ) ;
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	private:
		GenomePositionRange const m_range ;
	} ;
}

#endif
