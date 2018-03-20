
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_SNPTEST_STOPPING_CONDITION_HPP
#define METRO_SNPTEST_STOPPING_CONDITION_HPP

#include "appcontext/OstreamTee.hpp"

namespace metro {
	template< typename LogLikelihood >
	struct Snptest25StoppingCondition {
		Snptest25StoppingCondition(
			LogLikelihood const& ll,
			double tolerance,
			std::size_t max_iterations,
			appcontext::OstreamTee* log
		):
			m_ll( ll ),
			m_target_ll( -std::numeric_limits< double >::infinity() ),
			m_one_norm_of_first_derivative( std::numeric_limits< double >::infinity() ),
			m_current_ll( -std::numeric_limits< double >::infinity() ),
			m_tolerance( tolerance ),
			m_iteration( -1 ),
			m_max_iterations( max_iterations ),
			m_log( log )
		{}

		void set_target_loglikelihood( double const target ) {
			m_target_ll = target ;
		}

		bool operator()(
			LogLikelihood const& ll,
			Eigen::VectorXd const& parameters
		) {
			return operator()( parameters, ll.get_value_of_first_derivative() ) ;
		}

		bool operator()(
			Eigen::VectorXd const& parameters,
			Eigen::VectorXd const& first_derivative
		) {
			++m_iteration ;

			m_one_norm_of_first_derivative = first_derivative.array().abs().maxCoeff() ;

			double const previous_ll = m_current_ll ;
			m_current_ll = m_ll.get_value_of_function() ;

			if( m_current_ll > m_target_ll ) {
				m_target_ll = m_current_ll ;
			}

			// double one_norm = parameters.array().abs().maxCoeff() ; 
			Eigen::VectorXd difference ;
			if( m_current_params.size() > 0 ) {
				difference = ( parameters - m_current_params ).array().abs() ;
			} else {
				difference = LogLikelihood::Vector::Constant( parameters.size(), std::numeric_limits< double >::infinity() ) ;
			}

			//double const param_diff = ( m_current_params.size() == 0 ) ? 10000 : ( parameters - m_current_params ).array().abs().maxCoeff() ;
			m_current_params = parameters ;
			bool const converged = (
				std::abs( m_current_ll - previous_ll ) < m_tolerance        // loglikelihood has not changed much
				&& ( m_current_ll > ( m_target_ll - m_tolerance ) )         // loglikelihood at least as large as previously seen maximum
				&& ( m_one_norm_of_first_derivative < m_tolerance )			// derivative close to zero
			//	&& param_diff < m_tolerance									// parameters have not changed much.
			) ;

			if( m_log ) {
				(*m_log) << "Iteration: " << m_iteration << "; tolerance = " << m_tolerance << ":\n"
					<< "params: " << m_current_params.transpose() << "\n"
					<< "difference: " << difference.transpose() << " (rel diff  = " << ( difference.array() / ( parameters.array().abs() + 0.1 ) ).transpose() << ")\n"
					<< "current ll: " << m_current_ll << ", difference: " << std::abs( m_current_ll - previous_ll )
						<< ", rel. diff = " << std::abs( ( m_current_ll - previous_ll ) / m_current_ll ) << "\n"
					<< "target ll: " << m_target_ll << "\n"
					<< "first derivative: " << first_derivative.transpose() << "\n"
					//<< "second derivative:\n" << m_ll.get_value_of_second_derivative() << "\n"
					<< "converged: " << converged << ".\n" ;
			}
			return diverged() || converged ; 
		}

		bool diverged() const { return m_iteration >= m_max_iterations || m_current_ll != m_current_ll ; }
		std::size_t number_of_iterations() const { return m_iteration ; }

		private:
			LogLikelihood const& m_ll ;
			double m_target_ll ;
			double m_one_norm_of_first_derivative ;
			double m_current_ll ;
			Eigen::VectorXd m_current_params ;
			double const m_tolerance ;
			int m_iteration ;
			bool m_converged ;
			std::size_t const m_max_iterations ;
			appcontext::OstreamTee* m_log ;
	} ;
}

#endif
