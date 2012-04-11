
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef SNPIDMATCHES_CONDITION_HPP
#define SNPIDMATCHES_CONDITION_HPP

#include <set>
#include "GenRow.hpp"
#include "RowCondition.hpp"
#include "string_to_value_map.hpp"

struct SNPIDMatchesCondition: public RowCondition
{
public:
	SNPIDMatchesCondition( std::string const& expression ) ;
	bool check_if_satisfied( string_to_value_map const& ) const ;
	void format_to_stream( std::ostream& oStream ) const ;

private:
	void setup() ;

private:		
	std::string m_expression, m_prefix, m_suffix ;
	char m_wildcard_char ;
	
} ;



#endif

