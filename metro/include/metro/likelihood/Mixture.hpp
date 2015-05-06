
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
			
			void set_data( Matrix const& data ) {
				m_data = &data ;
				for( std::size_t i = 0; i < m_components.size(); ++i ) {
					m_components[i].set_data( data ) ;
				}
			}
			
			// Add a component with given name and weight to the mixture.
			// This takes ownership of the distribution
			template< typename DistributionPtr >
			void add_component(
				std::string const& name,
				Scalar const weight,
				DistributionPtr distribution
			) {
				Eigen::VectorXd new_parameters( m_parameters.size() + distribution->get_parameters().size() ) ;
				new_parameters.segment( 0, m_parameters.size() ) = m_parameters ;
				new_parameters.segment( m_parameters.size(), distribution->get_parameters().size() ) = distribution->get_parameters() ;

				m_component_names.push_back( name ) ;
				m_components.push_back( distribution ) ;
				m_weights.push_back( weight ) ;
				m_total_weight += weight ;
				m_parameters = new_parameters ;
			}
			
			void evaluate_at( Vector const& parameters ) {
				evaluate_at( parameters, DataSubset( DataRange( 0, m_data->rows() ))) ;
			}

			void evaluate_at( Vector const& parameters, DataSubset const& data_subset ) {
				m_parameters = parameters ;
				{
					int param_i = 0 ;
					for( std::size_t i = 0; i < m_components.size(); ++i ) {
						m_components[i].evaluate_at(
							parameters.segment(
								param_i,
								m_components[i].get_parameters().size()
							),
							data_subset
						) ;
						param_i += m_components[i].get_parameters().size() ;
					}
					assert( param_i == parameters.size() ) ;
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
				m_parameters = parameters ;
			}

			double get_value_of_function() const {
				return m_terms.sum() ;
			}
			
			void get_terms_of_function( MatrixRef result ) const {
				result = m_terms ;
			}

			Vector get_parameters() const {
				return m_parameters ;
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
				Vector m_parameters ;
				std::vector< Scalar > m_rescaled_weights ;
				Matrix m_component_terms ;
				Vector m_terms ;
				
		} ;
	}
}

#endif
