
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <string>
#include <cassert>
#include "genfile/GenomePosition.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"
#include "genfile/SNPPositionInRangeTest.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	SNPPositionInRangeTest::SNPPositionInRangeTest( GenomePositionRange const range ):
		m_range( range )
	{
	}
		
	bool SNPPositionInRangeTest::operator()( VariantIdentifyingData const& data ) const {
		return m_range.contains( data.get_position() ) ;
	}
	
	std::string SNPPositionInRangeTest::display() const {
		return "position in " + string_utils::to_string( m_range ) ;
	}
}
