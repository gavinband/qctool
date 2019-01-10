
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_SNPTEST_STOPPING_CONDITION_HPP
#define METRO_SNPTEST_STOPPING_CONDITION_HPP

//#include "appcontext/OstreamTee.hpp"

namespace metro {
	
	struct Trace {
	public:
		int iteration ;
		int max_iterations ;
		double tolerance ;

		double ll ;
		double target_ll ;
		double ll_difference ;
		double ll_rel_difference ;

		Eigen::VectorXd parameters ;
		Eigen::VectorXd difference ;
		Eigen::VectorXd first_derivative ;
		Eigen::MatrixXd second_derivative ;
		bool converged ;
		bool diverged ;
		
		std::ostream& operator<<( std::ostream& out ) {
			out << "Iteration: " << iteration << "; tolerance = " << tolerance << ":\n"
				<< "params: " << parameters.transpose() << "\n"
				<< "difference: " << difference.transpose() << " (rel diff  = " << ( difference.array() / ( parameters.array().abs() + 0.1 ) ).transpose() << ")\n"
				<< "current ll: " << ll << ", difference: " << ll_difference
				<< ", rel. diff = " << ll_rel_difference << "\n"
				<< "target ll: " << target_ll << "\n"
				<< "first derivative: " << first_derivative.transpose() << "\n"
				//<< "second derivative:\n" << m_ll.get_value_of_second_derivative() << "\n"
				<< "converged: " << converged << "\n"
				<< "diverged: " << diverged << ".\n" ;
			return out ;
		}
	} ;
	
	template< typename LogLikelihood >
	struct Snptest25StoppingCondition {
		Snptest25StoppingCondition(
			LogLikelihood const& ll,
			double tolerance,
			std::size_t max_iterations,
			Trace* trace
		):
			m_ll( ll ),
			m_target_ll( -std::numeric_limits< double >::infinity() ),
			m_one_norm_of_first_derivative( std::numeric_limits< double >::infinity() ),
			m_current_ll( -std::numeric_limits< double >::infinity() ),
			m_tolerance( tolerance ),
			m_iteration( -1 ),
			m_max_iterations( max_iterations ),
			m_trace( trace )
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

			if( m_trace ) {
				m_trace->iteration = m_iteration ;
				m_trace->max_iterations = m_max_iterations ;
				m_trace->parameters = m_current_params ;
				m_trace->ll = m_current_ll ;
				m_trace->target_ll = m_target_ll ;
				m_trace->ll_difference = std::abs( m_current_ll - previous_ll ) ;
				m_trace->ll_rel_difference = std::abs( ( m_current_ll - previous_ll ) / m_current_ll ) ;
				m_trace->first_derivative = first_derivative ;
				m_trace->second_derivative = m_ll.get_value_of_second_derivative() ;
				m_trace->converged = converged ;
				m_trace->diverged = diverged() ;
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
			Trace* m_trace ;
	} ;
}

#endif
