
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_IMPL_FILTER_GRAMMAR_HPP
#define GENFILE_IMPL_FILTER_GRAMMAR_HPP

#include <string>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/container.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SampleFilter.hpp"

namespace bf = boost::fusion ;

namespace genfile {
	namespace impl {
		SampleFilter::UniquePtr construct_filter( bf::vector< std::string, std::string, genfile::VariantEntry > const& value ) {
			assert(0) ;
		}

		SampleFilter::UniquePtr construct_between_filter(
			bf::vector< std::string, genfile::VariantEntry, genfile::VariantEntry > const& data
		) {
			SampleFilter::UniquePtr result(
				new VariableInInclusiveRangeSampleFilter( bf::at_c<0>( data ), bf::at_c<1>( data ), bf::at_c<2>( data ))
			) ;
			return result ;
		}

		using namespace boost::spirit::qi ;
		
		// This parser is adapted from a SQL grammaer by Andy Elvey
		// At http://boost-spirit.com/repository/applications/spirit_sql.zip
		template< typename Iterator >
		struct filter_parser : boost::spirit::qi::grammar< Iterator, SampleFilter::UniquePtr >
		{
			typedef std::string OP ;

			rule< Iterator, SampleFilter::UniquePtr > start ;
			rule< Iterator, SampleFilter::UniquePtr > op_clause ;
			rule< Iterator, SampleFilter::UniquePtr > between_clause ;
			rule< Iterator, SampleFilter::UniquePtr > in_clause ;
			rule< Iterator, SampleFilter::UniquePtr > NA_clause ;
			rule< Iterator, std::vector< genfile::VariantEntry > > value_list ;
			rule< Iterator, genfile::VariantEntry > value ;
			rule< Iterator, std::string > string_literal ;
			rule< Iterator, double > number ;
			symbols< char, std::string > op ;

    		filter_parser() : filter_parser::base_type( start, "filter_grammar" )
    		{
				char_		STAR('*');
				char_		COMMA(',');
				char_		LPAREN('(');
				char_		RPAREN(')'); 
				char_		SEMI(';');
				char_		LT('<');
				char_		GT('>');
				string		LE("<=");
				string		GE(">=");
				char_		EQUAL('=');
				string		NEQUAL("!=");
				
				typedef inhibit_case<strlit<> > token_t;

				token_t IN_			= as_lower["in"];
				token_t AND			= as_lower["and"];
				token_t OR			= as_lower["or"];
				token_t AS			= as_lower["as"];
				token_t IS			= as_lower["is"];
				token_t NOT			= as_lower["not"]; 
				token_t NA			= as_lower["NA"];  

				identifier =
					nocase
					[
						lexeme
						[
							+graph_p
						]
					]
				;

		        
				main_clause
					= longest[ var_op_clause | var_in_clause | var_between_clause | var_null_clause ] ;

				// Result of op_clause is
				// fusion_vector< a, b, c >
				// where a is std::vector< char >
				// b is op
				// c is genfile::VariantEntry.
			    op_clause
					=  varname >> op >> value[ _val = construct_filter( _1 ) ] ; 
				in_clause
					= varname >> !(NOT) >> IN_ >> value_list ; 
				between_clause
					= varname >> !(NOT) >> BETWEEN >> value >> AND >> value ; 
				NA_clause
					= varname >> IS >> !(NOT) >> NA ; 
				value_list %= LPAREN >> list_p( value, COMMA ) >> RPAREN ; 
				value %= string_literal | number ;
				string_literal %=
					lexeme
					[
					'\'' >> +( char_ - '\'' ) >> '\''
					]
					|
					lexeme
					[
					  '\"' >> +( char_ - '\"' ) >> '\"'
					]
				;
				number %= real_p ;
				op.add( "<", "<" )( "<=", "<=" )( ">", ">" )( ">=", ">=" )( "=", "=" )( "!=", "!=" ) ;
			} ;
		} ;
	}
}

#endif
