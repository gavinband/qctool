
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_rBF_HPP
#define METRO_rBF_HPP

#include <Eigen/Core>

namespace metro {
	// compute natural logarithm of rBF, the retrospective Bayes factor from 
	// Stephens and Balding, Nature Genetics 2009.
	// counts should be 2xn
	double compute_rBF( Eigen::MatrixXd const& counts, double const lambda ) ;
}

#endif
