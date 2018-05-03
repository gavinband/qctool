
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "Eigen/LU"
#include "metro/regression/Design.hpp"
#include "metro/regression/LogLikelihood.hpp"
#include "metro/regression/IndependentNormalWeightedLogLikelihood.hpp"

// #define DEBUG_NORMALWEIGHTEDLOGLIKELIHOOD 1

namespace metro {
	namespace regression {
		IndependentNormalWeightedLogLikelihood::UniquePtr IndependentNormalWeightedLogLikelihood::create(
			LogLikelihood::UniquePtr ll,
			std::vector< int > parameter_indices,
			std::vector< double > means,
			std::vector< double > variances
		) {
			return UniquePtr( new IndependentNormalWeightedLogLikelihood( ll, parameter_indices, means, variances ) ) ;
		}

		IndependentNormalWeightedLogLikelihood::IndependentNormalWeightedLogLikelihood(
			LogLikelihood::UniquePtr ll,
			std::vector< int > parameter_indices,
			std::vector< double > means,
			std::vector< double > variances
		):
			m_ll( ll ),
			m_parameter_indices( parameter_indices ),
			m_means( means ),
			m_variances( variances ),
			m_normalisation( eZeroAtMean ),
			m_constant( 0.0 ) // ll when at mean
		{
			if( m_normalisation == eZeroAtMean ) {
				m_constant = 0.0 ;
			} else if( m_normalisation == ePDF ) {
				assert(0) ; // This may not be correct yet
				m_constant = 0.0 ;
				for( std::size_t i = 0; i < m_parameter_indices.size(); ++i ) {
					m_constant -= std::log( std::sqrt( 2.0 * 3.14159265358979323846 * m_variances[i] )) ;
				}
			} else {
				assert(0) ;
			}
		}

		std::string IndependentNormalWeightedLogLikelihood::get_parameter_name( std::size_t i ) const {
			return m_ll->get_parameter_name(i) ;
		}
		IndependentNormalWeightedLogLikelihood::IntegerMatrix IndependentNormalWeightedLogLikelihood::identify_parameters() const {
			return m_ll->identify_parameters() ;
		}
		
		int IndependentNormalWeightedLogLikelihood::number_of_outcomes() const {
			return m_ll->number_of_outcomes() ;
		}
		
		void IndependentNormalWeightedLogLikelihood::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			m_ll->evaluate_at( parameters, numberOfDerivatives ) ;
			m_value_of_function = 0.0 ;
			m_value_of_first_derivative.setZero( parameters.size() ) ;
			m_value_of_second_derivative.setZero( parameters.size(), parameters.size() ) ;
			for( std::size_t i = 0; i < m_parameter_indices.size(); ++i ) {
				int const& parameter_index = m_parameter_indices[i] ;
				assert( parameter_index >= 0 && parameter_index < parameters.size() ) ;
				double const centralised = (parameters(parameter_index)-m_means[i]) ;

				m_value_of_function += -0.5 * (centralised * centralised) / m_variances[i] ;
				m_value_of_first_derivative( parameter_index ) = -centralised / m_variances[i] ;
				m_value_of_second_derivative( parameter_index, parameter_index )
					= -1.0 / m_variances[i] ;
			}
			m_value_of_function += m_constant ;
		}
	
		std::string IndependentNormalWeightedLogLikelihood::get_summary() const {
			std::ostringstream ostr ;
			ostr << "independently-normal-weighted:" << m_ll->get_summary() << "\n" ;
			ostr << "independently-normal-weighted: using the following prior parameters:\n" ;
			std::size_t max_label_length = 0 ;
			int const numberOfParameters = m_ll->identify_parameters().rows() ;
			for( int i = 0; i < numberOfParameters; ++i ) {
				max_label_length = std::max( max_label_length, m_ll->get_parameter_name(i).size() ) ;
			}
			ostr
				<< std::setw(3) << "" << "  "
				<< std::setw( max_label_length + 2 ) << ""
				<< ":"
				<< " " << std::setw(5) << "mean"
				<< " " << std::setw(8) << "variance"
				<< "\n" ;
			for( std::size_t i = 0; i < m_parameter_indices.size(); ++i ) {
				int const parameter_index = m_parameter_indices[i] ;
				ostr << std::setw(3) << parameter_index
					<< ": " << std::setw( max_label_length + 2 ) << m_ll->get_parameter_name(parameter_index) << ":" ;
					ostr
						<< " " << std::setw(5) << m_means[i]
						<< " " << std::setw(8) << m_variances[i]
						<< "\n" ;
			}
			return ostr.str() ;
		}
	}
}
