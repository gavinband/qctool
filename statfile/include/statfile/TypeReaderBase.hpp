#ifndef STATFILE_TYPEREADERBASE_HPP
#define STATFILE_TYPEREADERBASE_HPP

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
	struct TypeReaderBase: public TypeReaderBase< T2, T3, T4, T5, T6, T7, T8, T9, T10 >
	{
		using TypeReaderBase< T2, T3, T4, T5, T6, T7, T8, T9, T10 >::read_value ;
		virtual void read_value( T1& ) {} ;
	} ;

	template<>
	struct TypeReaderBase< empty, empty, empty, empty, empty, empty, empty, empty, empty, empty >
	{
		void read_value( empty& ) {} ;
		virtual ~TypeReaderBase() {} ;
	} ;
}

#endif
