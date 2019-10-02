
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <boost/noncopyable.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "Eigen/Core"
#include "Eigen/LU"
#include "metro/regression/Design.hpp"
#include "metro/regression/LogLikelihood.hpp"
#include "metro/regression/IndependentLogFWeightedLogLikelihood.hpp"

// #define DEBUG_NORMALWEIGHTEDLOGLIKELIHOOD 1

namespace metro {
	namespace regression {
		/*
		logF density
		this can be computed in R as follows:
		
		logistic <- function(x) { exp(x) / (1 + exp(x) )}
		dlogf <- function( x, nu1 = 1, nu2 = 1, normalisation = "zeroatmean" ) {
			alpha = nu1/2
			beta = nu2/2
			l = logistic((x))
			constant = 0
			if( normalisation == "zeroatmean" ) {
				constant = -(alpha+beta) * log(0.5)
			} else if( normalisation == "density" ) {
				constant = -lbeta( alpha, beta )
			} else {
				constant = 0
			}
			constant + l^alpha * ( 1-l)^ beta
		}
		dlogf.v2 <- function( x, nu1 = 1, nu2 = 1 ) {
			dbeta( logistic(x), shape1 = nu1/2, shape2 = nu2/2 )
		}
		*/
		IndependentLogFWeightedLogLikelihood::UniquePtr IndependentLogFWeightedLogLikelihood::create(
			LogLikelihood::UniquePtr ll,
			std::vector< int > parameter_indices,
			std::vector< double > nu1,
			std::vector< double > nu2
		) {
			return UniquePtr( new IndependentLogFWeightedLogLikelihood( ll, parameter_indices, nu1, nu2 ) ) ;
		}

		IndependentLogFWeightedLogLikelihood::IndependentLogFWeightedLogLikelihood(
			LogLikelihood::UniquePtr ll,
			std::vector< int > parameter_indices,
			std::vector< double > nu1,
			std::vector< double > nu2
		):
			m_ll( ll ),
			m_parameter_indices( parameter_indices ),
			m_alpha( nu1 ),
			m_beta( nu2 ),
			m_normalisation( eZeroAtMean )
		{
			assert( m_alpha.size() == m_parameter_indices.size() ) ;
			assert( m_beta.size() == m_parameter_indices.size() ) ;
			for( std::size_t i = 0; i < m_alpha.size(); ++i ) {
				m_alpha[i] = m_alpha[i]/2.0 ;
				m_beta[i] = m_beta[i]/2.0 ;
			}
			m_constant = 0.0 ;
			if( m_normalisation == eZeroAtMean ) {
				// Normalise so the log-likelihood value at zero is zero.
				// ll computes as log(0.5)*(alpha+beta) for each parameter.
				m_constant = 0.0  ;
				for( std::size_t i = 0; i < m_alpha.size(); ++i ) {
					m_constant -= (m_alpha[i] + m_beta[i]) * std::log(0.5) ;
				}
			} else if( m_normalisation == ePDF ) {
				assert(0) ; // this might not be correct yet
				using boost::math::lgamma ;
				for( std::size_t i = 0; i < m_alpha.size(); ++i ) {
					m_constant -= lgamma(m_alpha[i]) + lgamma(m_beta[i]) - lgamma(m_alpha[i] + m_beta[i] ) ;
				}
			} else {
				assert(0) ;
			}
		}

		std::string IndependentLogFWeightedLogLikelihood::get_parameter_name( std::size_t i ) const {
			return m_ll->get_parameter_name(i) ;
		}
		IndependentLogFWeightedLogLikelihood::IntegerMatrix IndependentLogFWeightedLogLikelihood::identify_parameters() const {
			return m_ll->identify_parameters() ;
		}
		
		int IndependentLogFWeightedLogLikelihood::number_of_parameters() const {
			return m_ll->number_of_parameters() ;
		}

		int IndependentLogFWeightedLogLikelihood::number_of_outcomes() const {
			return m_ll->number_of_outcomes() ;
		}
		
		void IndependentLogFWeightedLogLikelihood::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			m_ll->evaluate_at( parameters, numberOfDerivatives ) ;
			evaluate_impl( numberOfDerivatives ) ;
		}

		void IndependentLogFWeightedLogLikelihood::evaluate( int const numberOfDerivatives ) {
			m_ll->evaluate( numberOfDerivatives ) ;
			evaluate_impl( numberOfDerivatives ) ;
		}
	
		void IndependentLogFWeightedLogLikelihood::evaluate_impl( int const numberOfDerivatives ) {
			Vector const& parameters = m_ll->parameters() ;
			m_value_of_function = m_constant ;
			m_value_of_first_derivative.setZero( parameters.size() ) ;
			m_value_of_second_derivative.setZero( parameters.size(), parameters.size() ) ;

			for( std::size_t i = 0; i < m_parameter_indices.size(); ++i ) {
				int const& parameter_index = m_parameter_indices[i] ;
				assert( parameter_index >= 0 && parameter_index < parameters.size() ) ;
				double const x = parameters( parameter_index ) ;
				// logistic function evaluated at this value.
				double const l = 1.0 / (1.0 + std::exp(-x)) ;
				m_value_of_function += ( m_alpha[i] * std::log(l) + m_beta[i] * std::log( 1-l ) );
				m_value_of_first_derivative( parameter_index ) = m_alpha[i] - (m_alpha[i] + m_beta[i]) * l ;
				m_value_of_second_derivative( parameter_index, parameter_index ) = -(m_alpha[i] + m_beta[i]) * l * (1-l) ;
			}
		}

		IndependentLogFWeightedLogLikelihood::Vector IndependentLogFWeightedLogLikelihood::get_prior_mode() const {
			Vector result = Vector::Zero( m_value_of_first_derivative.size() ) ;
			for( std::size_t i = 0; i < m_parameter_indices.size(); ++i ) {
				// mode is mode of beta distribution with shpae1=1+alpha, shape2
				result(i) = m_alpha[i] / ( m_alpha[i] + m_beta[i] ) ;
			}
			return result ;
		}
	
		std::string IndependentLogFWeightedLogLikelihood::get_summary() const {
			std::ostringstream ostr ;
			ostr << m_ll->get_summary() << "\n" ;
			ostr << "  with priors:\n" ;
			std::size_t max_label_length = 0 ;
			int const numberOfParameters = m_ll->identify_parameters().rows() ;
			for( int i = 0; i < numberOfParameters; ++i ) {
				max_label_length = std::max( max_label_length, m_ll->get_parameter_name(i).size() ) ;
			}

			for( std::size_t i = 0; i < m_parameter_indices.size(); ++i ) {
				int const parameter_index = m_parameter_indices[i] ;
				ostr
					<< "  "
					<< std::setw( max_label_length)
					<< m_ll->get_parameter_name(parameter_index)
					<< " ~ logF( "
					<< (m_alpha[i]*2) << ", " << (m_beta[i]*2)
					<< " ).\n" ;
			}
#if 0			
			
			ostr
				<< std::setw(3) << "" << "  "
				<< std::setw( max_label_length + 2 ) << ""
				<< ":"
				<< " " << std::setw(5) << "alpha"
				<< " " << std::setw(5) << "beta"
				<< "\n" ;
			for( std::size_t i = 0; i < m_parameter_indices.size(); ++i ) {
				int const parameter_index = m_parameter_indices[i] ;
				ostr << std::setw(3) << parameter_index
					<< ": " << std::setw( max_label_length + 2 ) << m_ll->get_parameter_name(parameter_index) << ":" ;
					ostr
						<< " " << std::setw(5) << m_alpha[i]
						<< " " << std::setw(5) << m_beta[i]
						<< "\n" ;
			}
#endif
			return ostr.str() ;
		}
	}
}
