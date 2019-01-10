
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_REGRESSION_LOG_POSTERIOR_DENSITY_HPP
#define SNPTEST_REGRESSION_LOG_POSTERIOR_DENSITY_HPP

#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "metro/regression/Design.hpp"

#include "metro/regression/LogLikelihood.hpp"

namespace metro {
	namespace regression {
		/*
		*
		* DEPRECATED: use LogUnnormalisedPosterior instead
		*
		* Base class for classes implementing regression log posterior densities
		* These are prior / log-likelihood combinations that can unpack the contribution from likelihood and prior
		*/
		struct LogPosteriorDensity: public LogLikelihood
		{
			typedef std::auto_ptr< LogPosteriorDensity > UniquePtr ;
			virtual Vector get_prior_mode() const = 0 ;
			virtual Matrix get_loglikelihood_second_derivative() const = 0 ;
		} ;
	}
}

#endif
