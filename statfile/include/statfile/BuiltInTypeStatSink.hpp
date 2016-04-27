
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_BUILTINTYPESTATSINK_HPP
#define STATFILE_BUILTINTYPESTATSINK_HPP

#include <unistd.h>
#include <memory>
#include <string>
#include "genfile/Chromosome.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/StatSink.hpp"

namespace statfile {
	namespace {
		typedef StatSink< genfile::MissingValue, long unsigned int, int32_t, uint32_t, int64_t, uint64_t, std::string, double, genfile::Chromosome, genfile::GenomePosition > Base ;
	}
	struct BuiltInTypeStatSink: public Base
	{
		typedef std::auto_ptr< BuiltInTypeStatSink > UniquePtr ;
		static UniquePtr open( std::string const& filename ) ;
		using Base::write_value ;
		
		virtual void write_value( int32_t const& ) = 0 ;
		virtual void write_value( uint32_t const& ) = 0 ;
		virtual void write_value( int64_t const& ) = 0 ;
		virtual void write_value( uint64_t const& ) = 0 ;
		virtual void write_value( std::string const& ) = 0 ;
		virtual void write_value( double const& ) = 0 ;
		virtual void write_value( genfile::MissingValue const& ) = 0 ;
		virtual void write_value( genfile::Chromosome const& ) = 0 ;
		virtual void write_value( genfile::GenomePosition const& ) = 0 ;
	} ;
	
	struct NullBuiltInTypeStatSink: public ColumnNamingStatSink< BuiltInTypeStatSink >
	{
		static std::auto_ptr< BuiltInTypeStatSink > open() ;
		
		operator bool() const { return 0 ; }
		
		void write_value( int32_t const& ) { assert(0) ; }
		void write_value( uint32_t const& ) { assert(0) ; }
		void write_value( int64_t const& ) { assert(0) ; }
		void write_value( uint64_t const& ) { assert(0) ; }
		void write_value( std::string const& ) { assert(0) ; }
		void write_value( double const& ) { assert(0) ; }
		void write_value( genfile::MissingValue const& ) { assert(0) ; }
		void write_value( genfile::Chromosome const& ) { assert(0) ; }
		void write_value( genfile::GenomePosition const& ) { assert(0) ; }
	} ;
	
	struct TrivialBuiltInTypeStatSink: public ColumnNamingStatSink< BuiltInTypeStatSink >
	{
		void write_value( genfile::MissingValue const& ) {}
		void write_value( int32_t const& ) {}
		void write_value( uint32_t const& ) {}
		void write_value( int64_t const& ) {}
		void write_value( uint64_t const& ) {}
		void write_value( std::string const& ) {}
		void write_value( double const& ) {}
		void write_value( genfile::Chromosome const& ) {}
		void write_value( genfile::GenomePosition const& ) {}
		operator bool() const { return reinterpret_cast< void* >( const_cast< TrivialBuiltInTypeStatSink* >( this )) ; } ;
	} ;
}

#endif
