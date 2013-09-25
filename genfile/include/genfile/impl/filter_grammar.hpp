
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_IMPL_FILTER_GRAMMAR_HPP
#define GENFILE_IMPL_FILTER_GRAMMAR_HPP

#include <string>
#include "../../config.hpp"
#if HAVE_BOOST_SPIRIT
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_list.hpp>
#include <boost/spirit/include/qi_optional.hpp>
#include <boost/fusion/container.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SampleFilter.hpp"

namespace fusion = boost::fusion ;

namespace genfile {
	namespace impl {
		SampleFilter::UniquePtr construct_filter( fusion::vector< std::string, std::string, genfile::VariantEntry > const& data ) {
			assert(0) ;
		}

		SampleFilter::UniquePtr construct_in_filter( fusion::vector< std::string, std::vector< genfile::VariantEntry > > const& data ) {
			VariableInSetSampleFilter::UniquePtr result(
				new VariableInSetSampleFilter( fusion::at_c<0>( data ))
			) ;

			std::vector< genfile::VariantEntry > const& values = fusion::at_c<1>( data ) ;
			for( std::size_t i = 0; i < values.size(); ++i ) {
				result->add_level( values[i] ) ;
			}

			return SampleFilter::UniquePtr( result.release() ) ;
		}

		SampleFilter::UniquePtr construct_between_filter(
			fusion::vector< std::string, genfile::VariantEntry, genfile::VariantEntry > const& data
		) {
			SampleFilter::UniquePtr result(
				new VariableInInclusiveRangeSampleFilter( fusion::at_c<0>( data ), fusion::at_c<1>( data ), fusion::at_c<2>( data ))
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
			rule< Iterator, genfile::VariantEntry > identifier ;
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
				literal_char< char, true >		STAR('*') ;
				literal_char< char, true >		COMMA(',') ;
				literal_char< char, true >		LPAREN('(') ;
				literal_char< char, true >		RPAREN(')') ; 
				literal_char< char, true >		SEMI(';') ;
				literal_char< char, true >		LT('<') ;
				literal_char< char, true >		GT('>') ;
				literal_string< char, true >	LE("<=") ;
				literal_string< char, true >	GE(">=") ;
				literal_char< char, true >		EQUAL('=') ;
				literal_string< char, true >	NEQUAL("!=") ;

				literal_string< char, true >	IN_("IN") ;
				literal_string< char, true >	AND_("AND") ;
				literal_string< char, true >	OR_("OR") ;
				literal_string< char, true >	IS_("IS") ;
				literal_string< char, true >	NOT_("NOT") ;
				literal_string< char, true >	NA_("NA") ;
				
				start = op_clause | in_clause | between_clause | NA_clause ;

				// Result of op_clause is
				// fusion_vector< a, b, c >
				// where a is std::vector< char >
				// b is op
				// c is genfile::VariantEntry.
				op.add( "<", "<" )( "<=", "<=" )( ">", ">" )( ">=", ">=" )( "=", "=" )( "!=", "!=" ) ;
				identifier = +graph ;
			    op_clause			= ( identifier >> op >> value )[ _val = construct_filter( _1 ) ] ; 
				in_clause			= ( identifier >> IN_ >> value_list )[ _val = construct_in_filter( _1 ) ] ; 
				between_clause		= ( identifier >> IN_ >> '[' >> value >> ',' >> value >> ']' )[ _val = construct_between_filter( _1 ) ] ;
				NA_clause 			= identifier >> IS_ >> NA_ ;
				value_list			%= LPAREN >> ( value % ',' ) >> RPAREN ; 
				value				= string_literal | number ;
				string_literal		=
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
				number				= double_ ;
			} ;
		} ;
	}
}
#endif
#endif
