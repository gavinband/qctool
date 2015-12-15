
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <string>
#include <sstream>
#include <boost/bind.hpp>
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPIDInListTest.hpp"

namespace genfile {
	SNPIDInListTest::SNPIDInListTest( std::set< std::string > id_fields )
	: m_id_fields( id_fields )
	{}

	namespace {
		void insert_value( std::vector< std::string >* target, std::string const& id ) {
			target->push_back( id ) ;
		}
	}
	bool SNPIDInListTest::operator()( VariantIdentifyingData const& data ) const {
		std::vector< std::string > ids ;
		data.get_alleles( boost::bind( &insert_value, &ids, _1 )) ;
		for( std::size_t i = 0; i < ids.size(); ++i ) {
			if( m_id_fields.find( ids[i] ) != m_id_fields.end() ) {
				return true ;
			}
		}
		return false ;
	}
	
	std::string SNPIDInListTest::display() const {
		std::ostringstream ostr ;
		ostr << "SNPID in { " ;
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
