#ifndef STATFILE_BUILTINTYPESTATSINK_HPP
#define STATFILE_BUILTINTYPESTATSINK_HPP

#include <unistd.h>
#include <memory>
#include <string>
#include "genfile/Chromosome.hpp"
#include "statfile/StatSink.hpp"

namespace statfile {
	struct BuiltInTypeStatSink: public StatSink< int32_t, uint32_t, std::string, double, genfile::Chromosome >
	{
		typedef std::auto_ptr< BuiltInTypeStatSink > UniquePtr ;
		static UniquePtr open( std::string const& filename ) ;
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
