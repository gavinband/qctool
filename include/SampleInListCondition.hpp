
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef __GTOOL_SAMPLEINLIST_HPP__
#define __GTOOL_SAMPLEINLIST_HPP__

#include <set>
#include "GenRow.hpp"
#include "RowCondition.hpp"
#include "FileSet.hpp"
#include "string_to_value_map.hpp"

struct SampleInListCondition: public RowCondition
{
	SampleInListCondition( std::string filename ) ;
	bool check_if_satisfied( string_to_value_map const& ) const ;
	void format_to_stream( std::ostream& oStream ) const ;
	protected:
		FromFileSet< std::set< std::string > > m_id_list ;
		std::string m_filename ;
} ;



#endif

