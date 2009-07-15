
#ifndef __GTOOL_SNPINLIST_HPP__
#define __GTOOL_SNPINLIST_HPP__

#include <set>
#include "GenRow.hpp"
#include "RowCondition.hpp"
#include "FileSet.hpp"
#include "string_to_value_map.hpp"

struct SNPInListCondition: public RowCondition
{
	SNPInListCondition( std::string filename ) ;
	bool check_if_satisfied( string_to_value_map const& ) const ;
	void format_to_stream( std::ostream& oStream ) const ;
	protected:
		FromFileSet< std::set< std::string > > m_id_list ;
		std::string m_filename ;
} ;



#endif

