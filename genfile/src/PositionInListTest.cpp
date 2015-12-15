
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <string>
#include <sstream>
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/PositionInListTest.hpp"

namespace genfile {
	PositionInListTest::PositionInListTest( std::set< GenomePosition > positions )
		: m_positions( positions )
	{}

	bool PositionInListTest::operator()( VariantIdentifyingData const& data ) const {
		return m_positions.find( data.get_position() ) != m_positions.end() ;
	}

	std::string PositionInListTest::display() const {
		std::ostringstream ostr ;
		ostr << "Position in { " ;
		if( m_positions.size() <= 10 ) {
			for( std::set< GenomePosition >::const_iterator i = m_positions.begin(); i != m_positions.end(); ++i ) {
				if( i != m_positions.begin() ) {
					ostr << ", " ;
				}
				ostr << *i ;
			}
		}
		else {
			ostr << "set of " << m_positions.size() ;
		}
		ostr << " }" ;
		return ostr.str() ;
		
	}
}
