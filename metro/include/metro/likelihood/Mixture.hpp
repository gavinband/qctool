
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_MIXTURE_HPP
#define METRO_LIKELIHOOD_MIXTURE_HPP

#include <cmath>
#include <boost/ptr_container/ptr_vector.hpp>
#include <Eigen/Core>
#include "metro/IndependentObservationLogLikelihood.hpp"
#include "metro/DataRange.hpp"
#include "metro/DataSubset.hpp"
#include "metro/log_sum_exp.hpp"

// #define DEBUG_METRO_MIXTURE 1

namespace metro {
	namespace likelihood {
		template< typename Scalar, typename Vector, typename Matrix >
		struct Mixture: public metro::IndependentObservationLogLikelihood< Scalar, Vector, Matrix > {
			typedef metro::IndependentObservationLogLikelihood< Scalar, Vector, Matrix > Component ;
			typedef typename metro::IndependentObservationLogLikelihood< Scalar, Vector, Matrix >::MatrixRef MatrixRef ;
			
			Mixture( Matrix const& data ):
				m_data( &data ),
				m_total_weight(0)
			{}
			
			// Set data to a given matrix.
			// This will invalidate the parameters iff it does for any of the components.
			void set_data( Matrix const& data ) {
				m_data = &data ;
				for( std::size_t i = 0; i < m_components.size(); ++i ) {
					m_components[i].set_data( data ) ;
				}
			}
			
			// Add a component with given name and weight to the mixture.
			// This takes ownership of the distribution, which will be destroyed
			// when this class goes out of scope.
			template< typename DistributionPtr >
			void add_component(
				std::string const& name,
				Scalar const weight,
				DistributionPtr distribution
			) {
				m_component_names.push_back( name ) ;
				m_components.push_back( distribution ) ;
				m_weights.push_back( weight ) ;
				m_total_weight += weight ;
			}
			
			// Evaluate at a given set of parameters expressed as a single vector.
			// Parameters are taken in order of the components.
			void evaluate_at( Vector const& parameters ) {
				evaluate_at( parameters, DataSubset( DataRange( 0, m_data->rows() ))) ;
			}

			// Evaluate at a given set of parameters expressed as a single vector,
			// on a given subset of the data. Parameters are taken in order of the components.
			void evaluate_at( Vector const& parameters, DataSubset const& data_subset ) {
				{
					int param_i = 0 ;
					double total_weight = 0 ;
					for( std::size_t i = 0; i < m_components.size(); ++i ) {
						int const componentParameterSize = m_components[i].parameters().size() ;
						m_components[i].evaluate_at(
							parameters.segment(
								param_i,
								componentParameterSize
							),
							data_subset
						) ;
						if( i == ( m_components.size() - 1 ) ) {
							m_weights[i] = 1 - total_weight ;
							m_total_weight = 1 ;
						} else {
							m_weights[i] = parameters( param_i + componentParameterSize ) ;
							total_weight += m_weights[i] ;
						}
						param_i += m_components[i].parameters().size() + 1 ;
					}
					// sanity check
					assert( param_i == parameters.size() + 1 ) ;
				}

				m_component_terms.resize( m_data->rows(), m_components.size() ) ;
				m_component_terms.setZero() ;
				for( std::size_t i = 0; i < m_components.size(); ++i ) {
					m_components[i].get_terms_of_function( m_component_terms.col(i) ) ;
					// Re-weight according to cluster weights.
					m_component_terms.col(i) += Vector::Constant( m_component_terms.rows(), std::log( m_weights[i] / m_total_weight ) ) ;
#if DEBUG_METRO_MIXTURE
					std::cerr << "m_weights[" << i << "] = " << m_weights[i] / m_total_weight << "\n" ;
#endif
				}

#if DEBUG_METRO_MIXTURE
				std::cerr << "m_component_terms = \n"
					<< m_component_terms.block( 0, 0, 10, m_component_terms.cols() ) << "\n" ;
					std::cerr << "m_component_terms = \n"
						<< m_component_terms.block( 0, 0, 10, m_component_terms.cols() ) << "\n" ;
#endif
					

				// Compute it
				m_terms.setZero( m_component_terms.rows() ) ;
				rowwise_log_sum_exp( m_component_terms, &m_terms ) ;
			}

			double get_value_of_function() const {
				return m_terms.sum() ;
			}
			
			void get_terms_of_function( MatrixRef result ) const {
				result = m_terms ;
			}

			// Return the parameter vector.
			// This comes in the order:
			// Parameters for first component, (normalised) weight for first component
			// Parameters for second component, (normalised) weight for second component
			// ...
			// Parameters for last component.  No weight is included for last component as weights sum to one.
			Vector parameters() const {
				return compute_parameters() ;
			}

			std::string get_spec() const {
				return "Mixture" ;
			}

			Vector get_value_of_first_derivative() const {
				assert(0) ;
			}

			Matrix get_value_of_second_derivative() const {
				assert(0) ;
			}

		private:
			Matrix const* m_data ;
			boost::ptr_vector< Component > m_components ;
			std::vector< std::string > m_component_names ;
			std::vector< Scalar > m_weights ;
			Scalar m_total_weight ;
			std::vector< Scalar > m_rescaled_weights ;
			Matrix m_component_terms ;
			Vector m_terms ;

		private:
			Vector compute_parameters() const {
				if( m_components.size() == 0 ) {
					return Vector() ;
				} else {
					int numberOfParameters = 0 ;
					for( std::size_t i = 0; i < m_components.size(); ++i ) {
						numberOfParameters += m_components[i].parameters().size() + 1 ;
					}
					Eigen::VectorXd parameters( numberOfParameters - 1 ) ;
					int parameterIndex = 0 ;
					for( std::size_t i = 0; i < m_components.size(); ++i ) {
						Vector const& componentParams = m_components[i].parameters() ;
						parameters.segment( parameterIndex, componentParams.size() ) = componentParams ;
						if( i != ( m_components.size() - 1 ) ) {
							parameters( parameterIndex + componentParams.size() ) = m_weights[i] / m_total_weight ;
						}
						parameterIndex += componentParams.size() + 1 ;
					}
					return parameters ;
				}
			}
		} ;
	}
}

#endif
