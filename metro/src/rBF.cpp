
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/gamma.hpp>
#include "metro/rBF.hpp"

namespace metro {
	namespace {
		// compute the log of the density of the Dirichlet-multinomial distribution.
		// See http://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution
		template< typename Row1, typename Row2 >
		double log_dirichlet_multinomial( Row1 const& counts, Row2 const& lambdas ) {
			using boost::math::lgamma ;
			double result = 0.0 ;
			if( lambdas.maxCoeff() > 0 ) {
				result
					+= lgamma( lambdas.sum() )
					- lgamma( ( lambdas + counts ).sum() )
				;
			}

			for( int i = 0; i < counts.size(); ++i ) {
				if( lambdas(i) > 0 ) {
					result += lgamma( lambdas(i) + counts(i) ) - lgamma( lambdas(i) ) ;
				}
			}
			return result ;
		}
	}

	double compute_rBF( Eigen::MatrixXd const& counts, double const lambda ) {
		double const n = counts.sum() ;
		Eigen::RowVectorXd const colSums = counts.colwise().sum() ;
		Eigen::RowVectorXd const lambdas = lambda * colSums / n ;
		return std::exp(
			log_dirichlet_multinomial( counts.row(0), lambdas )
			+ log_dirichlet_multinomial( counts.row(1), lambdas )
			- log_dirichlet_multinomial( colSums, lambdas )
		) ;
	}
}

