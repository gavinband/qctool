#include <set>
#include <cassert>
#include <fstream>
#include <iostream>

#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SNPInListCondition.hpp"
#include "GenRowStatistics.hpp"
#include "string_utils.hpp"

SNPInListCondition::SNPInListCondition( std::string const& filename )
 : m_filenames( std::size_t(1), filename )
{
	setup() ;
}

SNPInListCondition::SNPInListCondition( std::vector< std::string > const& filenames )
 : m_filenames( filenames )
{
	setup() ;
}

void SNPInListCondition::setup() {
	for( std::size_t i = 0; i < m_filenames.size(); ++i ) {
		FromFileSet< std::set< std::string > > strings_in_file( m_filenames[i] ) ;
		m_id_list.insert( strings_in_file.begin(), strings_in_file.end() ) ;
	}
}

bool SNPInListCondition::check_if_satisfied( string_to_value_map const& statistics ) const {
	GenRowStatistics const* row_statistics_ptr = dynamic_cast< GenRowStatistics const* >( &statistics ) ;
	if( !row_statistics_ptr ) {
		throw ConditionException( "SNPInListCondition only supports GenRowStatistics." ) ;
	}
	return list_contains( row_statistics_ptr->row().SNPID() )
		|| list_contains( row_statistics_ptr->row().RSID() )
		|| list_contains( to_string( row_statistics_ptr->row().SNP_position() ) ) ;
}

bool SNPInListCondition::list_contains( std::string const& elt ) const {
	return m_id_list.find( elt ) != m_id_list.end() ;
}

void SNPInListCondition::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< "snp-in-list(" ;
	for( std::size_t i = 0; i < m_filenames.size(); ++i ) {
		if( i > 0 ) {
			oStream << ", " ;
		}
		oStream << m_filenames[i] ;
	}
	oStream << ")" ;
}