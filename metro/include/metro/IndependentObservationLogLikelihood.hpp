
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_INDEPENDENTOBSERVATION_LOGLIKELIHOOD_HPP
#define METRO_LIKELIHOOD_INDEPENDENTOBSERVATION_LOGLIKELIHOOD_HPP

#include <memory>
#include <boost/noncopyable.hpp>
#include <Eigen/Core>
#include "metro/LogLikelihood.hpp"

namespace metro {
	template< typename Scalar, typename Vector, typename Matrix >
	struct IndependentObservationLogLikelihood: public LogLikelihood< Scalar, Vector, Matrix > {
		typedef std::auto_ptr< IndependentObservationLogLikelihood > UniquePtr ;
		typedef typename Eigen::Ref< Matrix > MatrixRef ;

		// A function which returns the underlying terms of the log-likelihood.
		virtual void get_terms_of_function( MatrixRef result ) const = 0 ;
	} ;
}

#endif
