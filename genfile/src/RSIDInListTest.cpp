
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <string>
#include <sstream>
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/RSIDInListTest.hpp"

namespace genfile {
	RSIDInListTest::RSIDInListTest( std::set< std::string > id_fields )
		: m_id_fields( id_fields )
	{}

	bool RSIDInListTest::operator()( VariantIdentifyingData const& data ) const {
		return m_id_fields.find( data.get_rsid() ) != m_id_fields.end() ;
	}
	
	std::string RSIDInListTest::display() const {
		std::ostringstream ostr ;
		ostr << "RSID in { " ;
		if( m_id_fields.size() <= 10 ) {
			for( std::set< std::string >::const_iterator i = m_id_fields.begin(); i != m_id_fields.end(); ++i ) {
				if( i != m_id_fields.begin() ) {
					ostr << ", " ;
				}
				ostr << *i ;
			}
		}
		else {
			ostr << "set of " << m_id_fields.size() ;
		}
		ostr << " }" ;
		return ostr.str() ;
	}
}
