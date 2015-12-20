
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef POSITIONINLISTTEST_HPP
#define POSITIONINLISTTEST_HPP

#include <set>
#include <string>
#include "genfile/GenomePosition.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"

namespace genfile {
	struct PositionInListTest: public VariantIdentifyingDataTest
	{
		PositionInListTest( std::set< GenomePosition > id_fields ) ;
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	private:
		std::set< GenomePosition > const m_positions ;
	} ;
}

#endif
