
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <string>
#include <sstream>
#include "genfile/Chromosome.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/ChromosomeInSetTest.hpp"

namespace genfile {
	ChromosomeInSetTest::ChromosomeInSetTest( std::set< Chromosome > const& chromosomes ):
		m_chromosomes( chromosomes )
	{}
	
	bool ChromosomeInSetTest::operator()( VariantIdentifyingData const& data ) const {
		return m_chromosomes.find( data.get_position().chromosome() ) != m_chromosomes.end() ;
	}
	
	std::string ChromosomeInSetTest::display() const {
		std::ostringstream result ;
		result << "chromosome in { " ;
		std::size_t count = 0 ;
		for( std::set< Chromosome >::const_iterator i = m_chromosomes.begin(); i != m_chromosomes.end(); ++i ) {
			if( count++ > 0 ) {
				result << ", " ;
			}
			result << *i ;
		}
		result << " }" ;
		return result.str() ;
	}
}
