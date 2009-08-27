
#ifndef __GTOOL_SNPINLIST_HPP__
#define __GTOOL_SNPINLIST_HPP__

#include <set>
#include "GenRow.hpp"
#include "RowCondition.hpp"
#include "FileSet.hpp"
#include "string_to_value_map.hpp"

struct SNPInListCondition: public RowCondition
{
public:
	SNPInListCondition( std::string const& filename ) ;
	SNPInListCondition( std::vector< std::string > const& filenames ) ;
	bool check_if_satisfied( string_to_value_map const& ) const ;
	void format_to_stream( std::ostream& oStream ) const ;

private:
		void setup() ;
		bool list_contains( std::string const& ) const ;

private:		
		std::set< std::string > m_id_list ;
		std::vector< std::string > m_filenames ;
} ;



#endif

