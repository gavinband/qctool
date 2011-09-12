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
		typedef StatSink< genfile::MissingValue, long unsigned int, int32_t, uint32_t, int64_t, uint64_t, std::string, double, genfile::Chromosome > Base ;
	}
	struct BuiltInTypeStatSink: public Base
	{
		typedef std::auto_ptr< BuiltInTypeStatSink > UniquePtr ;
		static UniquePtr open( std::string const& filename ) ;
		
		using Base::write_value ;
		void write_value( genfile::VariantEntry const& value ) ;
	} ;
	
	struct NullBuiltInTypeStatSink: public ColumnNamingStatSink< BuiltInTypeStatSink >
	{
		static std::auto_ptr< BuiltInTypeStatSink > open() ;
		
		operator void*() const { return 0 ; }
	} ;
	
	struct TrivialBuiltInTypeStatSink: public ColumnNamingStatSink< BuiltInTypeStatSink >
	{
		operator void*() const { return reinterpret_cast< void* >( const_cast< TrivialBuiltInTypeStatSink* >( this )) ; } ;
	} ;
}

#endif
