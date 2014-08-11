
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_IMPL_FILTER_GRAMMAR_HPP
#define GENFILE_IMPL_FILTER_GRAMMAR_HPP

#include <string>
#include "../../config.hpp"
#if 0
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_list.hpp>
#include <boost/spirit/include/qi_optional.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/container.hpp>
#include <boost/variant/recursive_variant.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SampleFilter.hpp"

namespace fusion = boost::fusion ;

namespace genfile {
	namespace impl {
		struct sql_tree;

		typedef
			boost::variant<
				boost::recursive_wrapper<sql_tree>,
				genfile::VariantEntry
			>
		sql_tree_node ;

		struct sql_tree
		{
			std::string identifier ;
			std::string type ;
			std::string data ;
		} ;
	}
}

BOOST_FUSION_ADAPT_STRUCT(
	genfile::impl::sql_tree,
	( std::string, identifier )
	( std::string, type )
	( std::string, data )
)
// ;

namespace genfile {
	namespace impl {
		using namespace boost::spirit::qi ;
		
		// This parser is adapted from a SQL grammaer by Andy Elvey
		// At http://boost-spirit.com/repository/applications/spirit_sql.zip
		template< typename Iterator >
		struct filter_grammar : boost::spirit::qi::grammar< Iterator, sql_tree, ascii::space_type >
		{
			rule< Iterator, sql_tree, ascii::space_type > start ;
			symbols< char, std::string > op ;
/*
			rule< Iterator, sql_tree(), ascii::space_type > compound_clause ;
			rule< Iterator, sql_tree(), ascii::space_type > clause ;
			rule< Iterator, sql_tree(), ascii::space_type > op_clause ;
			rule< Iterator, sql_tree(), ascii::space_type > in_clause ;
*/
/*			
			rule< Iterator, std::string, ascii::space_type > identifier ;
			rule< Iterator, genfile::VariantEntry, ascii::space_type > value ;
			rule< Iterator, std::vector< genfile::VariantEntry >, ascii::space_type > value_list ;
			rule< Iterator, double, ascii::space_type > number ;
			rule< Iterator, std::string, ascii::space_type > string_literal ;
*/
			filter_grammar()
				: filter_grammar::base_type( start, "filter_grammar" )
			{
				using namespace boost::spirit::qi ;
				using boost::spirit::ascii::string ;
				
				op.add
					( "=", "=" )
					( "!=", "!=" )
					( "<", "<" )
					( ">", ">" )
					( "<=", "<=" )
					( ">=", ">=" )
				;

				start = string( "A" ) >> string( ">" ) >> string( "B" ) ;
//				op_clause			= ( identifier >> op >> value ) ; 									//	[ fusion::vector< string, string, boost::variant< string, double > > ]
//				identifier			= lexeme[ +graph ] ;												//	[ _val = _1 ] ;
//				value				= *(lexeme[ +graph ]) ;
				/*
				compound_clause 	= (
					clause >> *( lit( "AND" ) >> clause )
				) | (
					clause >> *( lit( "OR" ) >> clause )
				) ;
				clause = op_clause | in_clause ;

				in_clause			= ( identifier >> lit( "IN" ) >> value_list ) ;								//	[ fusion::vector< string, string, std::vector< boost::variant< string, double > > > ]

				// Result of op_clause is
				// fusion_vector< a, b, c >
				// where a is std::vector< char >
				// b is op
				// c is genfile::VariantEntry.
				op.add( "<", "<" )( "<=", "<=" )( ">", ">" )( ">=", ">=" )( "=", "=" )( "!=", "!=" ) ;
				
				value				= string_literal | number ;
				value_list			= lit( "(" ) >> value % ',' >> lit( ")" ) ;
				string_literal		=
					lexeme
					[
						'\'' >> ( +( char_ - '\'' )[ _val += _1 ] ) >> '\''
					]
					|
					lexeme
					[
						'\"' >> ( +( char_ - '\"' )[ _val += _1 ] ) >> '\"'
					]
				;
				number				= double_ ;
				*/
			} ;
		} ;
	}
}

#endif
#endif
