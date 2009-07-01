#include <set>
#include <cassert>
#include <fstream>

#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SampleInListCondition.hpp"
#include "SampleRowStatistics.hpp"

SampleInListCondition::SampleInListCondition( std::string filename )
 : m_id_list( filename ),
	m_filename( filename )
{}

bool SampleInListCondition::check_if_satisfied( GenotypeAssayStatistics const& statistics ) const {
	SampleRowStatistics const* row_statistics_ptr = dynamic_cast< SampleRowStatistics const* >( &statistics ) ;
	if( !row_statistics_ptr ) {
		throw ConditionException( "SampleInListCondition only supports SampleRowStatistics." ) ;
	}

	return ( m_id_list.find( row_statistics_ptr->row().ID1() ) != m_id_list.end() )
		|| ( m_id_list.find( row_statistics_ptr->row().ID2() ) != m_id_list.end() ) ;
}

void SampleInListCondition::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< "sample-in-list("
		<< m_filename
		<< ")" ;
}