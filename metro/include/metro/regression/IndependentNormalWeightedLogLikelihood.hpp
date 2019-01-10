
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_REGRESSION_INDEPENDENT_NORMAL_WEIGHTED_LOGLIKELIHOOD_HPP
#define SNPTEST_REGRESSION_INDEPENDENT_NORMAL_WEIGHTED_LOGLIKELIHOOD_HPP

#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "metro/regression/Design.hpp"
#include "metro/regression/LogLikelihood.hpp"
#include "metro/regression/LogPosteriorDensity.hpp"

namespace metro {
	namespace regression {
		// Represents a log-likelihood function weighted by independent
		// normal distributions on each parameter.
		// by default variance=infinity, i.e. no weighting, and weighting
		// must be set.
		struct IndependentNormalWeightedLogLikelihood: public LogPosteriorDensity {
		public:
			typedef std::auto_ptr< IndependentNormalWeightedLogLikelihood > UniquePtr ;
			enum Normalisation { ePDF, eZeroAtMean } ;
			static UniquePtr create(
				LogLikelihood::UniquePtr ll,
				std::vector< int > parameter_indices,
				std::vector< double > means,
				std::vector< double > variances
			) ;
		public:
			IndependentNormalWeightedLogLikelihood(
				LogLikelihood::UniquePtr ll,
				std::vector< int > parameter_indices,
				std::vector< double > means,
				std::vector< double > variances
			) ;
			
			regression::Design& design() const { return m_ll->design() ; }

			int number_of_parameters() const ;
			std::string get_parameter_name( std::size_t i ) const ;

			int number_of_outcomes() const ;
			// Return a lx2 matrix identifying the l parameters.
			// The row for each parameter contains the outcome level and design matrix column for that parameter.
			IntegerMatrix identify_parameters() const ;
			
			void evaluate_at( Point const& parameters, int const numberOfDerivatives = 2 ) ;

			Point const& get_parameters() const {
				return m_ll->get_parameters() ;
			}

			double get_value_of_function() const {
				return m_ll->get_value_of_function() + m_value_of_function ;
			}

			Vector get_value_of_first_derivative() const {
				return m_ll->get_value_of_first_derivative() + m_value_of_first_derivative ;
			}

			Matrix get_value_of_second_derivative() const {
				return m_ll->get_value_of_second_derivative() + m_value_of_second_derivative ;
			}
		
			Vector get_prior_mode() const ;

			Matrix get_prior_second_derivative() const {
				return m_value_of_second_derivative ;
			}

			Matrix get_loglikelihood_second_derivative() const {
				return m_ll->get_value_of_second_derivative() ;
			}
		
			std::string get_summary() const ;
			
		private:
			LogLikelihood::UniquePtr m_ll ;
			std::vector< int > const m_parameter_indices ;
			std::vector< double > const m_means ;
			std::vector< double > const m_variances ;
			Normalisation const m_normalisation ;
			double m_constant ;
			double m_value_of_function ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;
		} ;
	}
}
#endif
