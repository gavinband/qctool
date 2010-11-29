#ifndef NULLSTATSOURCE_HPP
#define NULLSTATSOURCE_HPP

#include <string>
#include <stddef.h>
#include "statfile/StatSource.hpp"

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
	struct NullStatSource: public StatSource< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >
	{
		std::size_t number_of_rows() const { return 0 ; }
		std::size_t number_of_columns() const { return 0 ; }
		std::vector< std::string > column_names() const {
			return std::vector< std::string >() ;
		}

	protected:

		void ignore_value() {}
		void ignore_all() {}
	} ;
	
}
#endif
