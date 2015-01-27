
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef NULLSTATSINK_HPP
#define NULLSTATSINK_HPP

#include "statfile/StatSink.hpp"

namespace statfile
{
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
	struct NullStatSink: public StatSink< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >
	{
		operator bool() const { return 0 ; }
		std::size_t number_of_columns() const { return 0 ; }
		std::vector< std::string > const& column_names() const {
			return m_column_names ;
		} ;

	private:
		void add_column_impl( std::string const& ) {} ;
		std::vector< std::string > m_column_names ;
	} ;
}

#endif
