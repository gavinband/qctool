
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_REGRESSION_COMBINED_LOGLIKELIHOOD_PRIOR_HPP
#define SNPTEST_REGRESSION_COMBINED_LOGLIKELIHOOD_PRIOR_HPP

#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "metro/regression/Design.hpp"
#include "metro/regression/LogLikelihood.hpp"
#include "metro/SmoothFunction.hpp"

namespace metro {
	namespace regression {
		/*
		* Base class for classes implementing regression log posterior densities
		* These are prior / log-likelihood combinations that can unpack the contribution from likelihood and prior
		*/
		struct LogUnnormalisedPosterior: public SmoothFunction
		{
		public:
			typedef std::auto_ptr< LogUnnormalisedPosterior > UniquePtr ;
			LogUnnormalisedPosterior( LogLikelihood& loglikelihood, SmoothFunction& posterior ) ;

		public:
			LogLikelihood const& ll() const ;
			SmoothFunction const& prior() const ;
			std::string get_summary() const ;
			
		public:
			int number_of_parameters() const { return m_ll.number_of_parameters() ; }
			
			void evaluate_at( Vector const& parameters, int const numberOfDerivatives = 2 ) ;
			void evaluate( int const numberOfDerivatives = 2 ) ;

			Scalar get_value_of_function() const {
				return m_ll.get_value_of_function() + m_prior.get_value_of_function() ;
			}
		
			Vector get_value_of_first_derivative() const {
				return m_ll.get_value_of_first_derivative() + m_prior.get_value_of_first_derivative() ;
			}
		
			Matrix get_value_of_second_derivative() const {
				return m_ll.get_value_of_second_derivative() + m_prior.get_value_of_second_derivative() ;
			}
		
			Vector const& parameters() const {
				return m_parameters ;
			}

			
		private:
			LogLikelihood& m_ll ;
			SmoothFunction& m_prior ;
			Vector m_parameters ;
		} ;
	}
}

#endif
