
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_PAIRWISE_CALL_COMPARER_HPP
#define QCTOOL_PAIRWISE_CALL_COMPARER_HPP

#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/VariantEntry.hpp"

struct PairwiseCallComparer {
public:
	typedef std::auto_ptr< PairwiseCallComparer > UniquePtr ;
	typedef boost::function< void ( std::string const&, genfile::VariantEntry const& ) > Callback ;
public:
	static UniquePtr create( std::string const& spec ) ;

	virtual ~PairwiseCallComparer() {}
	virtual void compare(
		Eigen::MatrixXd const& left,
		Eigen::MatrixXd const& right,
		Callback callback
	) const = 0 ;
} ;

#endif
