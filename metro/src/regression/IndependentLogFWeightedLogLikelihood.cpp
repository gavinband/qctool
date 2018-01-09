
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
#include "metro/regression::Design.hpp"
#include "metro/regression/LogLikelihood.hpp"
#include "metro/regression/IndependentLogFWeightedLogLikelihood.hpp"

// #define DEBUG_NORMALWEIGHTEDLOGLIKELIHOOD 1

namespace metro {
	namespace case_control {
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
			m_normalisation( eOneAtZero )
		{
			assert( m_alpha.size() == m_parameter_indices.size() ) ;
			assert( m_beta.size() == m_parameter_indices.size() ) ;
			for( std::size_t i = 0; i < m_alpha.size(); ++i ) {
				m_alpha[i] = m_alpha[i]/2.0 ;
				m_beta[i] = m_beta[i]/2.0 ;
			}
			m_constant = 0.0 ;
			if( m_normalisation == eOneAtZero ) {
				for( std::size_t i = 0; i < m_alpha.size(); ++i ) {
					m_constant = (m_alpha[i] + m_beta[i]) * std::log(0.5) ;
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
		
		int IndependentLogFWeightedLogLikelihood::number_of_outcomes() const {
			return m_ll->number_of_outcomes() ;
		}
		
		void IndependentLogFWeightedLogLikelihood::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			m_ll->evaluate_at( parameters, numberOfDerivatives ) ;
			m_value_of_function = 0.0 ;
			m_value_of_first_derivative.setZero( parameters.size() ) ;
			m_value_of_second_derivative.setZero( parameters.size(), parameters.size() ) ;

			for( std::size_t i = 0; i < m_parameter_indices.size(); ++i ) {
				int const& parameter_index = m_parameter_indices[i] ;
				assert( parameter_index >= 0 && parameter_index < parameters.size() ) ;
				double const x = parameters( parameter_index ) ;
				// logistic function evaluated at this value.
				double const l = 1.0 / (1.0 + std::exp(-x)) ;
				m_value_of_function += m_alpha[i] * std::log(l) + m_beta[i] * std::log( 1-l ) ;
				m_value_of_first_derivative( parameter_index ) = m_alpha[i] - (m_alpha[i] + m_beta[i]) * l ;
				m_value_of_second_derivative( parameter_index, parameter_index ) = -(m_alpha[i] + m_beta[i]) * l * (1-l) ;
			}
		}
	
		std::string IndependentLogFWeightedLogLikelihood::get_summary() const {
			std::ostringstream ostr ;
			ostr << "independently-logF-weighted:" << m_ll->get_summary() << "\n" ;
			ostr << "independently-logF-weighted: using the following prior weights:\n" ;
			std::size_t max_label_length = 0 ;
			int const numberOfParameters = m_ll->identify_parameters().rows() ;
			for( int i = 0; i < numberOfParameters; ++i ) {
				max_label_length = std::max( max_label_length, m_ll->get_parameter_name(i).size() ) ;
			}
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
			return ostr.str() ;
		}
	}
}
