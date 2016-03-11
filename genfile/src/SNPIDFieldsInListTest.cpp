
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <string>
#include "genfile/VariantIdentifyingDataTest.hpp"
#include "genfile/SNPIDFieldsInListTest.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	SNPIDFieldsInListTest::SNPIDFieldsInListTest( std::set< std::string > id_fields )
	: m_id_fields( id_fields )
	{}

	namespace {
		void insert_value( std::vector< std::string >* target, std::string const& id ) {
			target->push_back( id ) ;
		}
	}

	bool SNPIDFieldsInListTest::operator()( VariantIdentifyingData const& data ) const {
		if( m_id_fields.find( data.get_primary_id() ) != m_id_fields.end() ) {
			return true ;
		}
		std::vector< std::string > ids ;
		data.get_alleles( boost::bind( &insert_value, &ids, _1 )) ;
		for( std::size_t i = 0; i < ids.size(); ++i ) {
			if( m_id_fields.find( ids[i] ) != m_id_fields.end() ) {
				return true ;
			}
		}
		return false ;
	}
	
	std::string SNPIDFieldsInListTest::display() const {
		std::string result = "SNPID or RSID in { " ;
		if( m_id_fields.size() <= 10 ) {
			for( std::set< std::string >::const_iterator i = m_id_fields.begin(); i != m_id_fields.end(); ++i ) {
				if( i != m_id_fields.begin() ) {
					result += ", " ;
				}
				result += *i ;
			}
		}
		else {
			result += "set of " + string_utils::to_string( m_id_fields.size() ) ;
		}
		return result + " }";
	}
}
