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
	
	bool ChromosomeInSetTest::operator()(
		std::string,
		std::string,
		GenomePosition pos,
		std::string,
		std::string
	) const {
		return m_chromosomes.find( pos.chromosome() ) != m_chromosomes.end() ;
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
