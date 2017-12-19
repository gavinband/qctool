
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
#include "metro/case_control/RegressionDesign.hpp"
#include "metro/case_control/LogLikelihood.hpp"
#include "metro/case_control/NormalWeightedLogLikelihood.hpp"

// #define DEBUG_NORMALWEIGHTEDLOGLIKELIHOOD 1

namespace metro {
	namespace case_control {
		namespace impl {
			double compute_mvn_constant( int const k, double determinant ) {
				return -0.5 * k * std::log( 2 * 3.141592654 ) - 0.5 * std::log( determinant ) ;
			}
		}

		NormalWeightedLogLikelihood::UniquePtr NormalWeightedLogLikelihood::create(
			LogLikelihood::UniquePtr ll,
			Vector const mean,
			Matrix const covariance,
			Normalisation normalisation
		) {
			return UniquePtr( new NormalWeightedLogLikelihood( ll, mean, covariance, normalisation ) ) ;
		}

		NormalWeightedLogLikelihood::NormalWeightedLogLikelihood(
			LogLikelihood::UniquePtr ll,
			Vector const mean,
			Matrix const covariance,
			Normalisation normalisation
		):
			m_ll( ll ),
			m_mean( mean ),
			m_covariance( covariance ),
			m_normalisation( normalisation ),
			m_solver( covariance ),
			m_inverse_covariance( m_solver.solve( Matrix::Identity( mean.size(), mean.size() ))),
			m_constant( impl::compute_mvn_constant( m_mean.size(), m_covariance.determinant() ))
		{}

		void NormalWeightedLogLikelihood::set_predictor_levels(
			Matrix const& levels,
			Matrix const& probabilities,
			std::vector< metro::SampleRange > const& included_samples
		) { 
			m_ll->set_predictor_levels( levels, probabilities, included_samples ) ;
		}

		std::string NormalWeightedLogLikelihood::get_parameter_name( std::size_t i ) const {
			return m_ll->get_parameter_name(i) ;
		}
		NormalWeightedLogLikelihood::IntegerMatrix NormalWeightedLogLikelihood::identify_parameters() const {
			return m_ll->identify_parameters() ;
		}
		
		int NormalWeightedLogLikelihood::number_of_outcomes() const {
			return m_ll->number_of_outcomes() ;
		}
		
		void NormalWeightedLogLikelihood::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			m_ll->evaluate_at( parameters, numberOfDerivatives ) ;
			Vector solved = m_solver.solve( parameters - m_mean ) ;
			m_log_density = -0.5 * ((parameters - m_mean).transpose() * solved)(0) ;
			m_value_of_first_derivative = -solved ;
			m_value_of_second_derivative = -m_inverse_covariance ;

#if DEBUG_NORMALWEIGHTEDLOGLIKELIHOOD
			std::cerr << "parameters = " << parameters.transpose() << ".\n" ;
			std::cerr << "constant = " << m_constant << ".\n" ;
			std::cerr << "solved = \n" << solved << ".\n" ;
			std::cerr << "log density = " << m_log_density << "\n" ;
#endif
		}
	
		std::string NormalWeightedLogLikelihood::get_summary() const {
			std::ostringstream ostr ;
			ostr << "normal-weighted:" << m_ll->get_summary() << "\n" ;
			ostr << "normal-weighted: using the following prior covariance:\n" ;
			std::size_t max_label_length = 0 ;
			for( int i = 0; i < m_covariance.rows(); ++i ) {
				max_label_length = std::max( max_label_length, m_ll->get_parameter_name(i).size() ) ;
			}
			ostr
				<< std::setw(3) << "" << "  "
				<< std::setw( max_label_length + 2 ) << "parameter"
				<< ":" ;
			for( int j = 0; j < m_covariance.cols(); ++j ) {
				ostr << " " << std::setw(5) << (j+1) ;
			}
			ostr << "\n" ;
			for( int i = 0; i < m_covariance.rows(); ++i ) {
				ostr << std::setw(3) << i+1
					<< ": " << std::setw( max_label_length + 2 ) << m_ll->get_parameter_name(i) << ":" ;
				for( int j = 0; j < m_covariance.cols(); ++j ) {
					ostr << " " << std::setw(5) << std::setprecision(4) << m_covariance(i,j) ;
				}
				ostr << "\n" ;
			}
			return ostr.str() ;
		}
	}
}
