
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_TYPEWRITERBASE_HPP
#define STATFILE_TYPEWRITERBASE_HPP

#include "statfile/types.hpp"

namespace statfile {
	
	template<
		typename T1 = empty,
		typename T2 = empty,
		typename T3 = empty,
		typename T4 = empty,
		typename T5 = empty,
		typename T6 = empty,
		typename T7 = empty,
		typename T8 = empty,
		typename T9 = empty,
		typename T10 = empty
	>
	struct TypeWriterBase: public TypeWriterBase< T2, T3, T4, T5, T6, T7, T8, T9, T10 >
	{
		using TypeWriterBase< T2, T3, T4, T5, T6, T7, T8, T9, T10 >::write_value ;
		virtual void write_value( T1 const& ) {} ;
	} ;

	template<>
	struct TypeWriterBase< empty, empty, empty, empty, empty, empty, empty, empty, empty, empty >
	{
		virtual void write_value( empty const& ) {} ;
		virtual ~TypeWriterBase() {} ;
	} ;
}

#endif
