
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SAMPLEIDSINK_HPP
#define SAMPLEIDSINK_HPP

#include <vector>
#include <string>
#include "statfile/BuiltInTypeStatSink.hpp"
#include "ObjectSink.hpp"
#include "SampleRow.hpp"
#include "GToolException.hpp"
#include "string_utils/string_utils.hpp"

struct SampleIDSinkException: public std::exception
{
	char const* what() const throw() { return "SampleIDSinkException" ; }
} ;

// A SampleRow sink which just outputs sample ids, one per line,
// to the given stream.
struct SampleIDSink: public ObjectSink< SampleRow >
{
	SampleIDSink( std::string const& filename )
	: 	m_sink( statfile::BuiltInTypeStatSink::open( filename ) )
	{
		assert( m_sink.get() ) ;
		(*m_sink) | "ID_1" | "ID_2" ;
	}

	SampleIDSink& write( SampleRow const& row ) {
		(*m_sink) << row.ID1() << row.ID2() << statfile::end_row() ;
		return *this ;
	}

	operator bool() const {
		return (*m_sink) ;
	}

private:

	statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
} ;

#endif
