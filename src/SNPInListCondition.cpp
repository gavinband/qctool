#include <set>
#include <cassert>
#include <fstream>

#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SNPInListCondition.hpp"
#include "GenRowStatistics.hpp"

SNPInListCondition::SNPInListCondition( std::string filename )
 : m_id_list( filename ),
	m_filename( filename )
{}

bool SNPInListCondition::check_if_satisfied( GenotypeAssayStatistics const& statistics ) const {
	GenRowStatistics const* row_statistics_ptr = dynamic_cast< GenRowStatistics const* >( &statistics ) ;
	if( !row_statistics_ptr ) {
		throw ConditionException( "SNPInListCondition only supports GenRowStatistics." ) ;
	}
	return m_id_list.find( row_statistics_ptr->row().SNPID() ) != m_id_list.end() ;
}

void SNPInListCondition::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< "snp-in-list("
		<< m_filename
		<< ")" ;
}