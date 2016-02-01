
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <cassert>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"
#include "genfile/SNPInSetTest.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	SNPInSetTest::SNPInSetTest( std::set< VariantIdentifyingData > const& snps, VariantIdentifyingData::CompareFields const& comparer ):
	 	m_snps( snps.begin(), snps.end(), comparer )
	{}

	bool SNPInSetTest::operator()( VariantIdentifyingData const& data ) const {
		return( m_snps.find( data ) != m_snps.end() ) ;
	}
	
	std::string SNPInSetTest::display() const {
		std::ostringstream ostr ;
		ostr << "SNP (" + m_snps.key_comp().get_summary() + ") in { " ;
		int count = 0;
		for( std::set< VariantIdentifyingData >::const_iterator i = m_snps.begin(); i != m_snps.end() && count < 3; ++i, ++count ) {
			if( i != m_snps.begin() ) {
				ostr << ", " ;
			}
				ostr << (*i) ;
		}
		if( m_snps.size() > 3 ) {
			ostr << "...(+" << m_snps.size() << " others)" ;
		}
		ostr << " }" ;
		return ostr.str() ;
	}
}
