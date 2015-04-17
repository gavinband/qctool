
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
#include "metro/DataRange.hpp"

// #define DEBUG_MULTIVARIATE_T 1

namespace metro {
	namespace likelihood {
		template< typename Scalar, typename Vector, typename Matrix >
		struct MultivariateT: public metro::LogLikelihood< Scalar, Vector, Matrix > {
			MultivariateT( Matrix const& data, double const degrees_of_freedom ):
				m_pi( 3.141592653589793238462643383279502884 ),
				m_nu( degrees_of_freedom ),
				m_data( data ),
				m_data_range( std::vector< DataRange >( 1, DataRange( 0, m_data.rows() ))),
				// parameters are:
				// degrees of freedom
				// mean 
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

			MultivariateT( Matrix const& data, std::vector< DataRange > const& data_range, double const degrees_of_freedom ):
				m_pi( 3.141592653589793238462643383279502884 ),
				m_nu( degrees_of_freedom ),
				m_data( data ),
				m_data_range( data_range ),
				// parameters are:
				// degrees of freedom
				// mean 
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
					for( std::size_t i = 0; i < m_data_range.size(); ++i ) {
						DataRange const& range = m_data_range[i] ;
#if DEBUG_MULTIVARIATE_T
					std::cerr << "metro::likelihood::MultivariateT::get_value_of_function():\n"
						<< " Adding " << range << "\n" ;
#endif				
						result +=
							( range.size() * m_kappa )
							- 0.5 * ( m_Z2.segment( range.begin(), range.size() ).sum() + range.size() * m_log_determinant ) ;
					}
				} else {
					for( std::size_t i = 0; i < m_data_range.size(); ++i ) {
						DataRange const& range = m_data_range[i] ;
						result += -(( m_nu + m_p ) * (
							m_Z2.segment( range.begin(), range.size() )
							+ Vector::Constant( range.size(), m_nu )
						).array().log().sum() ) / 2.0 ;
						result += range.size() * ( m_kappa - 0.5 * m_log_determinant ) ;
					}
				}
				return result ;
			}

			Vector get_value_of_first_derivative() const {
				assert(0) ;
			}

			Matrix get_value_of_second_derivative() const {
				assert(0) ;
			}


			// Fit multivariate T by EM until the loglikelihood increases by less than the given amount.
			template< typename StoppingCondition >
			bool estimate_by_em(
				StoppingCondition& stopping_condition
			) {
				Matrix regularising_sigma = Matrix::Constant( m_p, m_p, 0 ) ;
				double regularising_weight = 0 ;
				return estimate_by_em(
					stopping_condition,
					regularising_sigma,
					regularising_weight
				) ;
			}

			// Fit multivariate T by EM until the stopping condition is satisfied.
			// Algorithm details are from Nadarajah & Kotz, "Estimation methods for the Multivariate t Distribution.", p.103
			// This function includes a regularising variance-covariance matrix, and weight
			// to prevent the fit from becoming degenerate.
			//
			// StoppingCondition must support operator() with result convertible to bool.
			// a value of true means stop, otherwise continue.
			// It must also be callable as stopping_condition.converged(), which this function uses
			// to return a value to the caller.
			template< typename StoppingCondition >
			bool estimate_by_em(
				StoppingCondition& stopping_condition,
				Matrix const& regularising_sigma,
				double regularising_weight
			) {
				assert( regularising_sigma.rows() == m_p ) ;
				assert( regularising_sigma.cols() == m_p ) ;
				assert( regularising_weight >= 0.0 ) ;

				// Start with unit weights, giving MVN estimate
				Vector weights = Vector::Constant( m_data.rows(), 1 ) ;
				Vector mean = compute_weighted_mean( weights, m_data_range ) ;
				Matrix sigma = compute_weighted_regularised_sigma( weights, mean, regularising_sigma, regularising_weight, m_data_range ) ;

				evaluate_at( mean, sigma ) ;
				double loglikelihood = get_value_of_function() ;
				
#if DEBUG_MULTIVARIATE_T
				std::cerr << "metro::likelihood::MultivariateT::estimate_by_em(): start: params = "
					<< get_parameters().transpose() << ", ll = " << get_value_of_function() << ".\n" ;
#endif				

				// If nu = ∞ we are at the MLE already, so bail out...
				if( m_nu == std::numeric_limits< double >::infinity() ) {
					// Multivariate normal.  Quit right now.
					return true ;
				}

				// ..otherwise let's EM it.
				std::size_t iteration = 0 ;
				while( !stopping_condition( loglikelihood ) ) {
					// compute weights
					// Vector of weights is given as
					// (nu+p) / nu + (x_i-mean)^t R^-1 ( x_i - mean ).
					// Our x_i - mean_i is stored in a single row of m_mean_centred_data.
					Vector weights = (
						m_mean_centred_data.array() * m_ldlt.solve( m_mean_centred_data.transpose() ).transpose().array()
					).rowwise().sum() ;
					weights += Vector::Constant( m_data.rows(), m_nu ) ;
					weights.array() = weights.array().inverse() * ( m_nu + m_p ) ;

					// compute new parameter estimates
					// these are
					// mean = sum( w_i x_i ) / sum( w_i )
					mean = compute_weighted_mean( weights, m_data_range ) ;
					sigma = compute_weighted_regularised_sigma( weights, mean, regularising_sigma, regularising_weight, m_data_range ) ;
					evaluate_at( mean, sigma ) ;
					loglikelihood = get_value_of_function() ;

#if DEBUG_MULTIVARIATE_T
					std::cerr << "metro::likelihood::MultivariateT::estimate_by_em(): after iteration "
						<< iteration << ": params = " << get_parameters().transpose()
						<< ", ll = " << loglikelihood << ", weights = " << weights.transpose() << ".\n" ;
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
			double const m_nu ;
			Matrix const& m_data ;
			std::vector< DataRange > m_data_range ;
			double const m_p ;
			double const m_kappa ;
			Vector m_parameters ;
			Matrix m_sigma ;
			Vector m_mean ;

			Eigen::LDLT< Matrix > m_ldlt ;
			double m_log_determinant ;
			Matrix m_mean_centred_data ;
			Vector m_Z2 ;
			
			Vector compute_weighted_mean(
				Vector const& weights,
				std::vector< DataRange > const& data_range
			) const {
				Vector result = Vector::Zero( m_p ) ;
				double total_weight = 0 ;
				for( std::size_t i = 0; i < data_range.size(); ++i ) {
					DataRange const& range = data_range[i] ;
					result += (
						(
							weights.segment( range.begin(), range.size() ).asDiagonal()
							* m_data.block( range.begin(), 0, range.size(), m_p )
						).colwise().sum()
					).transpose() ;

					total_weight += weights.segment( range.begin(), range.size() ).sum() ;
				}
				result /= total_weight ;
				return result ;
			}
			
			Matrix compute_weighted_regularised_sigma(
				Vector const& weights,
				Vector const& mean,
				Matrix const& regularising_sigma,
				double const regularising_weight,
				std::vector< DataRange > const& data_range
			) const {
				Matrix mean_centred_data = m_data.rowwise() - mean.transpose() ;
				// sigma = 1/N sum ( w_i (x_i-mean)(x_i-mean)^t)
				// we store x_i as a row not a column, so transposes go the opposite way.
				Matrix result = Matrix::Zero( m_p, m_p ) ;
				double n = 0 ;
				for( std::size_t i = 0; i < data_range.size(); ++i ) {
					DataRange const& range = data_range[i] ;
					result += (
						(
							mean_centred_data.block( range.begin(), 0, range.size(), m_p ).transpose()
							* weights.segment( range.begin(), range.size() ).asDiagonal()
							* mean_centred_data.block( range.begin(), 0, range.size(), m_p )
						)
					) ;
					n += range.size() ;
				}
				result += regularising_weight * regularising_sigma ;
				result /= ( n + regularising_weight ) ;
				return result ;
			}
		} ;
	}
}

#endif
