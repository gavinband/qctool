
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
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	class SampleFilter: public boost::noncopyable
	{
	public:
		typedef std::auto_ptr< SampleFilter > UniquePtr ;
		static UniquePtr create( std::string const& spec ) ;
		
		typedef Eigen::MatrixXd DetailMatrix ;
		typedef Eigen::Block< DetailMatrix > DetailBlock ;
	public:
		SampleFilter() ;
		virtual ~SampleFilter() ;

		virtual void summarise( std::ostream& ) const = 0 ;
		virtual std::size_t number_of_clauses() const { return 1 ; }
		virtual bool test(
			genfile::CohortIndividualSource const&,
			std::size_t i,
			DetailBlock* detail = 0
		) const = 0 ;
		void test(
			genfile::CohortIndividualSource const&,
			boost::function< void( std::size_t ) > callback,
			DetailMatrix* detail = 0
		) const ;
	} ;

	std::ostream& operator<< ( std::ostream& oStream, SampleFilter const& filter ) ;
}

#endif
