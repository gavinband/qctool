
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "test_case.hpp"
#include "genfile/impl/filter_grammar.hpp"
#include <cassert>

#if 1
AUTO_TEST_CASE( filter_grammar_test1 ) {
	
	using namespace boost::spirit::qi ;
	
	std::vector< std::string > sqls ;
	sqls.push_back( "A = 1" ) ;
	sqls.push_back( "A = \"string\"" ) ;
	sqls.push_back( "A = 'string'" ) ;
	sqls.push_back( "A != 1" ) ;
	sqls.push_back( "A != \"string\"" ) ;
	sqls.push_back( "A != 'string'" ) ;

	for( std::size_t i = 0; i < sqls.size(); ++i ) {
		genfile::impl::sql_tree tree ;
		genfile::impl::filter_grammar< std::string::const_iterator > grammar ;
		TEST_ASSERT( parse( sqls[i].begin(), sqls[i].end(), grammar, boost::spirit::ascii::space, tree ) ) ;
	}
}


#endif
