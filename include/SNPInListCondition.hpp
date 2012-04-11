
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNPINLIST_HPP
#define QCTOOL_SNPINLIST_HPP

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

	static std::string make_key( std::string const& SNPID, std::string const& RSID, uint32_t SNP_position ) ;

private:
		void setup() ;
		bool list_contains( std::string const& ) const ;

		bool file_appears_to_be_plain( std::string const& filename ) ;
		void load_from_gen_file( std::string const& filename ) ;
		void load_from_plain_file( std::string const& filename ) ;

private:		
		std::set< std::string > m_id_list ;
		std::vector< std::string > m_filenames ;
} ;



#endif

