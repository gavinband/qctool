
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_GENOTYPE_TABULATING_CALL_COMPARER_HPP
#define QCTOOL_GENOTYPE_TABULATING_CALL_COMPARER_HPP

#include <string>
#include <map>
#include <Eigen/Core>
#include <boost/function.hpp>
#include "genfile/VariantEntry.hpp"
#include "PairwiseCallComparer.hpp"

struct GenotypeTabulatingCallComparer: public PairwiseCallComparer {

	GenotypeTabulatingCallComparer() ;

	void compare(
		Eigen::MatrixXd const& left,
		Eigen::MatrixXd const& right,
		Callback
	) const ;
	
private:
	double const m_threshhold ;
} ;

#endif
