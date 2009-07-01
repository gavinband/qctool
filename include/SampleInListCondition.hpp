
#ifndef __GTOOL_SAMPLEINLIST_HPP__
#define __GTOOL_SAMPLEINLIST_HPP__

#include <set>
#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "RowCondition.hpp"
#include "FileUtil.hpp"

struct SampleInListCondition: public RowCondition
{
	SampleInListCondition( std::string filename ) ;
	bool check_if_satisfied( GenotypeAssayStatistics const& ) const ;
	void format_to_stream( std::ostream& oStream ) const ;
	protected:
		FromFileSet< std::set< std::string > > m_id_list ;
		std::string m_filename ;
} ;



#endif

