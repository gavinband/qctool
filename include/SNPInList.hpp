
#ifndef __GTOOL_SNPINLIST_HPP__
#define __GTOOL_SNPINLIST_HPP__

#include <set>
#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "RowCondition.hpp"
#include "FileUtil.hpp"

struct SNPInList: public RowCondition
{
	SNPInList( std::string filename ) ;
	bool check_if_satisfied( GenRow const& genRow, GenotypeAssayStatistics const * ) const ;

	protected:
		FromFileSet< std::set< std::string > > m_id_list ;
} ;



#endif

