
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPIDINLISTTEST_HPP
#define SNPIDINLISTTEST_HPP

#include <set>
#include <string>
#include "genfile/VariantIdentifyingDataTest.hpp"

namespace genfile {
	struct SNPIDInListTest: public VariantIdentifyingDataTest
	{
		SNPIDInListTest( std::set< std::string > id_fields ) ;
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	private:
		std::set< std::string > m_id_fields ;
	} ;
}

#endif
