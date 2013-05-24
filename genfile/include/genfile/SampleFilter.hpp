
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef GENFILE_SAMPLE_FILTER_HPP
#define GENFILE_SAMPLE_FILTER_HPP

#include <memory>
#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	class SampleFilter: public boost::noncopyable
	{
	public:
		typedef std::auto_ptr< SampleFilter > UniquePtr ;
		static UniquePtr create( std::string const& spec ) ;
	public:
		SampleFilter() ;
		virtual ~SampleFilter() ;

		virtual void summarise( std::ostream& ) const = 0 ;

		void compute_failed_samples( genfile::CohortIndividualSource const&, boost::function< void( std::size_t ) > callback ) const ;
		virtual bool test( genfile::CohortIndividualSource const&, std::size_t i ) const = 0 ;
	} ;

	std::ostream& operator<< ( std::ostream& oStream, SampleFilter const& filter ) ;
}
#endif