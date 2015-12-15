
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef RSIDINLISTTEST_HPP
#define RSIDINLISTTEST_HPP

#include <set>
#include <string>
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	struct RSIDInListTest: public SNPIdentifyingDataTest
	{
		RSIDInListTest( std::set< std::string > id_fields ) ;
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	private:
		std::set< std::string > m_id_fields ;
	} ;
}

#endif
