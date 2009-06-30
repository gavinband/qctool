#include <set>
#include <cassert>
#include <fstream>

#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SNPInListCondition.hpp"

SNPInListCondition::SNPInListCondition( std::string filename )
 : m_id_list( filename ),
	m_filename( filename )
{}

bool SNPInListCondition::check_if_satisfied( GenRow const& genRow, GenotypeAssayStatistics const* ) const {
	return m_id_list.find( genRow.SNPID() ) != m_id_list.end() ;
}

void SNPInListCondition::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< "snp-in-list("
		<< m_filename
		<< ")" ;
}