
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_PRODUCT_OF_MULTINOMIALS_HPP
#define METRO_LIKELIHOOD_PRODUCT_OF_MULTINOMIALS_HPP

#include <cmath>
#include "metro/LogLikelihood.hpp"
#include "metro/likelihood/Multinomial.hpp"

namespace metro {
	namespace likelihood {
		template< typename Scalar, typename Vector, typename Matrix >
		struct ProductOfMultinomials: public metro::LogLikelihood< Scalar, Vector, Matrix > {
			ProductOfMultinomials( Matrix const& data ):
				m_data( data )
			{
				setup( data, Matrix::Zero( data.rows(), data.cols() )) ;
			}
			
			ProductOfMultinomials( Matrix const& data, Matrix const& parameters ):
				m_data( data )
			{
				setup( data, parameters ) ;
			}

			void setup( Matrix const& data, Matrix const& parameters ) {
				assert( parameters.rows() == data.rows() && parameters.cols() == data.cols() ) ;
				for( int i = 0; i < data.rows(); ++i ) {
					m_multinomials.push_back( Multinomial< Scalar, Vector, Matrix >( data.row( i ), parameters.row( i ) )) ;
				}
			}

			void evaluate_at( Matrix const& parameters ) {
				assert( parameters.cols() == m_data.cols() && parameters.rows() == m_data.rows() ) ;
				for( std::size_t i = 0; i < m_multinomials.size(); ++i ) {
					m_multinomials[i].evaluate_at( parameters.row( i )) ;
				}
			}

			void evaluate_at( Vector const& parameters ) {
				assert( parameters.size() == m_data.rows() * m_data.cols() ) ;
				for( std::size_t i = 0; i < m_multinomials.size(); ++i ) {
					m_multinomials[i].evaluate_at( parameters.segment( i * m_data.cols(), m_data.cols() ) ) ;
				}
			}

			double get_value_of_function() const {
				double result = 0 ;
				for( std::size_t i = 0; i < m_multinomials.size(); ++i ) {
					result += m_multinomials[i].get_value_of_function() ;
				}
				return result ;
			}

			// For the purposes of derivatives we treat the parameters as a vector.
			// The ith block of N components are for the ith multinomial in the product.
			Vector get_value_of_first_derivative() const {
				Vector result = Vector::Zero( m_data.rows() * m_data.cols() ) ;
				for( std::size_t i = 0; i < m_multinomials.size(); ++i ) {
					result.segment( m_data.cols() * i, m_data.cols() ) = m_multinomials[i].get_value_of_first_derivative() ;
				}
				return result ;
			}

			// For the purposes of derivatives we treat the parameters as a vector.
			// The ith block of N components are for the ith multinomial in the product.
			Matrix get_value_of_second_derivative() const {
				Matrix result = Matrix::Zero( m_data.rows() * m_data.cols(), m_data.rows() * m_data.cols() ) ;
				for( std::size_t i = 0; i < m_multinomials.size(); ++i ) {
					result.block( m_data.cols() * i, m_data.cols() * i, m_data.cols(), m_data.cols() ) = m_multinomials[i].get_value_of_second_derivative() ;
				}
				return result ;
			}

			Vector get_MLE() const {
				Vector result( m_data.rows() * m_data.cols() ) ;
				for( std::size_t i = 0; i < m_multinomials.size(); ++i ) {
					result.segment( m_data.cols() * i, m_data.cols() ) = m_multinomials[i].get_MLE() ;
				}
				return result ;
			}

			std::string get_spec() const { return "ProductOfMultinomials" ; }
			
			Vector get_parameters() const {
				Vector result( m_data.rows() * m_data.cols() ) ;
				for( std::size_t i = 0; i < m_multinomials.size(); ++i ) {
					result.segment( m_data.cols() * i, m_data.cols() ) = m_multinomials[i].get_parameters() ;
				}
				return result ;
			}
		private:
			Matrix const m_data ;
			std::vector< Multinomial< Scalar, Vector, Matrix > > m_multinomials ;
		} ;
	}
}

#endif
