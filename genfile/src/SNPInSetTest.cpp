
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <cassert>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPInSetTest.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	SNPInSetTest::SNPInSetTest( std::set< SNPIdentifyingData > const& snps ):
	 	m_snps( snps )
	{}

 	bool SNPInSetTest::operator()( SNPIdentifyingData const& data ) const {
		return( m_snps.find( data ) != m_snps.end() ) ;
	}

	bool SNPInSetTest::operator()(
		std::string SNPID,
		std::string RSID,
		GenomePosition position,
		std::string alleleA,
		std::string alleleB
	) const {
		return operator()( 
			SNPIdentifyingData(
				SNPID, RSID, position, alleleA, alleleB
			)
		) ;
	}
	
	std::string SNPInSetTest::display() const {
		std::ostringstream ostr ;
		ostr << "SNP in { " ;
		if( m_snps.size() <= 3 ) {
			for( std::set< SNPIdentifyingData >::const_iterator i = m_snps.begin(); i != m_snps.end(); ++i ) {
				if( i != m_snps.begin() ) {
					ostr << ", " ;
				}
				ostr << (*i) ;
			}
		}
		else {
			ostr << "set of " + m_snps.size() ;
		}
		ostr << " }" ;
		return ostr.str() ;
	}
}
