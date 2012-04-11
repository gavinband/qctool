
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_BUILTINTYPESTATSOURCE_HPP
#define STATFILE_BUILTINTYPESTATSOURCE_HPP

#include <memory>
#include <string>
#include <vector>
#include "genfile/wildcard.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/Error.hpp"
#include "statfile/StatSource.hpp"

namespace statfile {
	struct BuiltInTypeStatSource: public StatSource< int32_t, uint32_t, std::string, double, genfile::Chromosome, char >
	{
		typedef StatSource< int32_t, uint32_t, std::string, double, genfile::Chromosome, char > Base ;
		typedef std::auto_ptr< BuiltInTypeStatSource > UniquePtr ;
		static UniquePtr open( std::string const& filename ) ;
		static UniquePtr open( std::vector< genfile::wildcard::FilenameMatch > const& filenames ) ;
		
	public:
		using Base::read_value ;
		// Convenience functions.
		// These act to expand the number of readable types.
		void read_value( genfile::Chromosome& chromosome ) ;
		void read_value( char& c ) ;
	} ;
	
	struct NullBuiltInTypeStatSource: public BuiltInTypeStatSource
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
