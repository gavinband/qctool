
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_LOGLIKELIHOOD_HPP
#define METRO_LIKELIHOOD_LOGLIKELIHOOD_HPP

#include <memory>
#include <boost/noncopyable.hpp>
#include <Eigen/Core>

namespace metro {
	template< typename Scalar, typename Vector, typename Matrix >
	struct LogLikelihood: public boost::noncopyable {
		typedef std::auto_ptr< LogLikelihood > UniquePtr ;
		virtual ~LogLikelihood() {}
		
		virtual Scalar get_value_of_function() const = 0 ;
		virtual Vector get_value_of_first_derivative() const = 0 ;
		virtual Matrix get_value_of_second_derivative() const = 0 ;
		
		virtual std::string get_spec() const = 0 ;
	} ;
}

#endif
