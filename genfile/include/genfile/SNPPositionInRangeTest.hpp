
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
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	// Test if the given SNP is in the given (closed) range of positions.
	struct SNPPositionInRangeTest: public SNPIdentifyingDataTest
	{
		SNPPositionInRangeTest( GenomePositionRange const range ) ;
		
		bool operator()(
			std::string,
			std::string,
			GenomePosition position,
			std::string,
			std::string
		) const ;

		std::string display() const ;
	private:
		GenomePositionRange const m_range ;
	} ;
}

#endif
