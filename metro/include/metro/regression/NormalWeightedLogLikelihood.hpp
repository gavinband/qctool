
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_REGRESSION_NORMAL_WEIGHTED_LOGLIKELIHOOD_HPP
#define SNPTEST_REGRESSION_NORMAL_WEIGHTED_LOGLIKELIHOOD_HPP

#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "metro/regression/Design.hpp"
#include "metro/regression/LogLikelihood.hpp"

namespace metro {
	namespace regression {
		// Represents a log-likelihood function weighted by independent
		// normal distributions on each parameter.
		// by default variance=infinity, i.e. no weighting, and weighting
		// must be set.
		struct NormalWeightedLogLikelihood: public LogLikelihood {
		public:
			typedef std::auto_ptr< NormalWeightedLogLikelihood > UniquePtr ;
			enum Normalisation { eNoConstantTerm = 0, eWithConstantTerm = 1 } ;
			static UniquePtr create(
				LogLikelihood::UniquePtr ll,
				Vector const mean,
				Matrix const covariance,
				Normalisation normalisation = eWithConstantTerm
			) ;
				
		public:
			NormalWeightedLogLikelihood(
				LogLikelihood::UniquePtr ll,
				Vector const mean,
				Matrix const covariance,
				Normalisation normalisation = eWithConstantTerm
			) ;
			
			regression::Design& design() const { return m_ll->design() ; }

			int number_of_parameters() const ;
			std::string get_parameter_name( std::size_t i ) const ;

			int number_of_outcomes() const ;

			// Return a lx2 matrix identifying the l parameters.
			// The row for each parameter contains the outcome level and design matrix column for that parameter.
			IntegerMatrix identify_parameters() const ;
			
			void evaluate_at( Point const& parameters, int const numberOfDerivatives = 2 ) ;
			void evaluate( int const numberOfDerivatives = 2 ) ;

			Point const& parameters() const {
				return m_ll->parameters() ;
			}

			double get_value_of_function() const {
				return m_ll->get_value_of_function() + m_log_density + ( (m_normalisation == eWithConstantTerm) ? m_constant : 0.0 ) ;
			}

			Vector get_value_of_first_derivative() const {
				return m_ll->get_value_of_first_derivative() + m_value_of_first_derivative ;
			}

			Matrix get_value_of_second_derivative() const {
				return m_ll->get_value_of_second_derivative() + m_value_of_second_derivative ;
			}
		
			std::string get_summary() const ;
				
		private:
			LogLikelihood::UniquePtr m_ll ;
			Vector const m_mean ;
			Matrix const m_covariance ;
			Normalisation const m_normalisation ;
			Eigen::LDLT< Matrix > m_solver ;
			Matrix const m_inverse_covariance ;
			double const m_constant ;
			double m_log_density ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;
			
		private:
			void evaluate_impl( int const numberOfDerivatives = 2 ) ;
		} ;
	}
}
#endif
