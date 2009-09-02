#include <set>
#include <cassert>
#include <fstream>

#include "SampleRow.hpp"
#include "SampleInListCondition.hpp"
#include "SampleRowStatistics.hpp"
#include "string_to_value_map.hpp"

SampleInListCondition::SampleInListCondition( std::string filename )
 : m_id_list( filename ),
	m_filename( filename )
{}

bool SampleInListCondition::check_if_satisfied( string_to_value_map const& statistics ) const {
	SampleRow const* row_ptr = dynamic_cast< SampleRow const* >( &statistics ) ;
	assert( row_ptr ) ;

	return ( m_id_list.find( row_ptr->ID1() ) != m_id_list.end() )
		|| ( m_id_list.find( row_ptr->ID2() ) != m_id_list.end() ) ;
}

void SampleInListCondition::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< "sample-in-list("
		<< m_filename
		<< ")" ;
}