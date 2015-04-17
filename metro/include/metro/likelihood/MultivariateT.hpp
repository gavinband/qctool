
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_MULTIVARIATE_T_HPP
#define METRO_LIKELIHOOD_MULTIVARIATE_T_HPP

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include "metro/LogLikelihood.hpp"

// #define DEBUG_MULTIVARIATE_T 1

namespace metro {
	namespace likelihood {
		template< typename Scalar, typename Vector, typename Matrix >
		struct MultivariateT: public metro::LogLikelihood< Scalar, Vector, Matrix > {
			MultivariateT( Matrix const& data, double const degrees_of_freedom ):
				m_pi( 3.141592653589793238462643383279502884 ),
				m_data( data ),
				// parameters are:
				// degrees of freedom
				// mean 
				m_nu( degrees_of_freedom ),
				m_p( data.cols() ),
				m_parameters( Vector::Zero( m_p + ( m_p * ( m_p + 1 ) / 2 ) )),
				// Compute contribution to log-likelihood that does not depend
				// on the parameters.
				m_kappa(
					degrees_of_freedom == std::numeric_limits< double >::infinity()
						? (
							-0.5 * m_p * std::log( 2 * m_pi )
						) : (
							lgamma( ( m_nu + m_p ) / 2.0 )
							- ( m_p / 2 ) * std::log( m_pi * m_nu )
							- lgamma( m_nu / 2.0 )
							+ (( m_nu + m_p ) * std::log( m_nu ) / 2.0 )
						)
				)
			{
			}

			MultivariateT( MultivariateT const& other ):
				m_data( other.m_data ),
				m_parameters( other.m_parameters )
			{}

			MultivariateT& operator=( MultivariateT const& other ) {
				m_data = other.m_data ;
				m_parameters = other.m_parameters ;
				return *this ;
			}

			// Evaluate at a packed set of parameters
			// These are:
			// The p values of the mean, followed by
			// The p (p+1)/2 entries of the lower diagonal of the sigma matrix (in column-major order.)
			void evaluate_at( Vector const& parameters ) {
				assert( parameters.size() == m_parameters.size() ) ;
				// Unpack parameters into mean and sigma
				m_mean = parameters.head( m_data.cols() ) ;
				int parameter_i = m_data.cols() ;
				for( int col = 0; col < m_data.cols(); ++col ) {
					m_sigma.col( col ).segment( col, m_sigma.cols() - col ) = parameters.segment( parameter_i, m_data.cols() - col ) ;
					parameter_i += m_data.cols() - col ;
				}
				evaluate_at( m_mean, m_sigma ) ;
			}

			void evaluate_at( Vector const& mean, Matrix const& sigma ) {
				assert( mean.size() == m_data.cols() ) ;
				assert( sigma.rows() == sigma.cols() ) ;
				assert( sigma.rows() == m_data.cols() ) ;

				m_mean = mean ;
				m_sigma = sigma ;

				// Pack parameters into the parameter vector
				m_parameters.segment( 0, m_mean.size() ) = m_mean ;
				int parameter_i = m_mean.size() ;
				for( int col = 0; col < m_data.cols(); ++col ) {
					m_parameters.segment( parameter_i, m_sigma.cols() - col ) = m_sigma.col( col ).segment( col, m_sigma.cols() - col ) ;
					parameter_i += m_data.cols() - col ;
				}
				
				// Uses only lower-diagonal.
				m_ldlt.compute( m_sigma ) ;
				m_log_determinant = m_ldlt.vectorD().array().log().sum() ;
				m_mean_centred_data = m_data.rowwise() - m_mean.transpose() ;
				m_Z2 = (
					m_mean_centred_data.array()
					* ( m_ldlt.solve( m_mean_centred_data.transpose() ).transpose().array() )
				).rowwise().sum() ;

#if DEBUG_MULTIVARIATE_T > 2
				std::cerr << "metro::likelihood::MultivariateT::evaluate_at():\n"
					<< " data = " << m_data << ".\n"
					<< " determinant( sigma ) = " << std::exp( m_log_determinant ) << "\n"
					<< " mean centred data = " << m_mean_centred_data.transpose() << ".\n" ;
#endif				
			}

			double get_value_of_function() const {
				double result = 0 ;
				if( m_nu == std::numeric_limits< double >::infinity() ) {
					// Multivariate normal.
					result = -0.5 * ( m_Z2.sum() + m_data.rows() * m_log_determinant ) + ( m_data.rows() * m_kappa ) ;
				} else {
					result = -(( m_nu + m_p ) * ( m_Z2 + Vector::Constant( m_data.rows(), m_nu ) ).array().log().sum() ) / 2.0 ;
					result -= m_data.rows() * 0.5 * m_log_determinant ;
					result += m_data.rows() * m_kappa ;
#if DEBUG_MULTIVARIATE_T > 2
					std::cerr << "metro::likelihood::MultivariateT::get_value_of_function():\n"
						<< " kappa = " << m_kappa << ".\n"
						<< " terms = " << (( m_nu * m_p ) * m_Z2.array().log().transpose() ) / 2.0 << ".\n" ;
#endif				
				}
				return result ;
			}

			Vector get_value_of_first_derivative() const {
				assert(0) ;
			}

			Matrix get_value_of_second_derivative() const {
				assert(0) ;
			}

			// Fit multivariate T by EM until the loglikelihood
			// increases by less than the given amount.
			// Algorithm details are from Nadarajah & Kotz, "Estimation methods for the Multivariate t Distribution.", p.103
			template< typename StoppingCondition >
			bool estimate_by_em( StoppingCondition& stopping_condition ) {
				// Start with the MLE multivariate normal estimate
				Vector mu = m_data.colwise().sum().transpose() / m_data.rows() ;
				Eigen::MatrixXd Z = m_data.rowwise() - mu.transpose() ;
				Matrix sum_of_squares = ( Z.transpose() * Z ) ;
				Matrix sigma = sum_of_squares / m_data.rows() ;
				evaluate_at( mu, sigma ) ;
				double loglikelihood = get_value_of_function() ;
				
#if DEBUG_MULTIVARIATE_T
				std::cerr << "metro::likelihood::MultivariateT::estimate_by_em(): start: params = "
					<< get_parameters().transpose() << ", ll = " << get_value_of_function() << ".\n" ;
#endif				

				if( m_nu == std::numeric_limits< double >::infinity() ) {
					// Multivariate normal.  Quit right now.
					return true ;
				}

				Vector weights = Vector::Constant( m_data.rows(), 1 ) ;
				std::size_t iteration = 0 ;
				while( !stopping_condition( loglikelihood ) ) {
					// compute weights
					// Vector of weights is given as
					// (nu+p) / nu + (x_i-mu)^t R^-1 ( x_i - mu ).
					// Our x_i - mu_i is stored in a single row of m_mean_centred_data.
				
					Vector weights = (
						m_mean_centred_data.array() * m_ldlt.solve( m_mean_centred_data.transpose() ).transpose().array()
					).rowwise().sum() ;
					weights += Vector::Constant( m_data.rows(), m_nu ) ;
					weights.array() = weights.array().inverse() * ( m_nu + m_p ) ;

					// compute new parameter estimates
					// these are
					// mu = sum( w_i x_i ) / sum( w_i )
					Eigen::VectorXd new_mu = (( weights.asDiagonal() * m_data ).colwise().sum() / weights.sum() ).transpose() ;
					m_mean_centred_data = m_data.rowwise() - new_mu.transpose() ;
					// sigma = 1/N sum ( w_i (x_i-mu)(x_i-mu)^t)
					// we store x_i as a row not a column, so transposes go the opposite way.
					sigma = ( m_mean_centred_data.transpose() * weights.asDiagonal() * m_mean_centred_data ) / m_data.rows() ;
					evaluate_at( new_mu, sigma ) ;
					loglikelihood = get_value_of_function() ;

#if DEBUG_MULTIVARIATE_T
					std::cerr << "metro::likelihood::MultivariateT::estimate_by_em(): after iteration "
						<< iteration << ": params = " << get_parameters().transpose() << ", ll = " << loglikelihood << ".\n" ;
#endif				
					++iteration ;
				}
				return stopping_condition.converged() ;
			}

			Vector const& get_parameters() const { return m_parameters ; }
			double get_degrees_of_freedom() const { return m_nu ; }
			Vector const& get_mean() const { return m_mean ; }
			Matrix const& get_sigma() const { return m_sigma ; }
			Matrix const& get_data() const { return m_data ; }
			
			std::string get_spec() const { return "MultivariateT" ; }

		private:
			double const m_pi ;
			Matrix const m_data ;
			double const m_nu ;
			double const m_p ;
			double const m_kappa ;
			Vector m_parameters ;
			Matrix m_sigma ;
			Vector m_mean ;

			Eigen::LDLT< Matrix > m_ldlt ;
			double m_log_determinant ;
			Eigen::MatrixXd m_mean_centred_data, m_Z2, m_Z3 ;
		} ;
	}
}

#endif
