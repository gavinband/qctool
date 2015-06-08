
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_MULTINOMIAL_HPP
#define METRO_LIKELIHOOD_MULTINOMIAL_HPP

#include <cmath>
#include "metro/LogLikelihood.hpp"

namespace metro {
	namespace likelihood {
		template< typename Scalar, typename Vector, typename Matrix >
		struct Multinomial: public metro::LogLikelihood< Scalar, Vector, Matrix > {
			Multinomial( Vector const& data ):
				m_counts( data ),
				m_parameters( Vector::Zero( data.size() ) )
			{
			}

			Multinomial( Vector const& counts, Vector const& parameters ):
				m_counts( counts ),
				m_parameters( parameters )
			{
				assert( parameters.size() == counts.size() ) ;
			}
			
			Multinomial( Multinomial const& other ):
				m_counts( other.m_counts ),
				m_parameters( other.m_parameters )
			{}

			Multinomial& operator=( Multinomial const& other ) {
				m_counts = other.m_counts ;
				m_parameters = other.m_parameters ;
				return *this ;
			}

			void evaluate_at( Vector const& parameters, DataSubset const& subset = DataRange( 0, 1 ) ) {
				assert( subset.number_of_subranges() == 1 && subset[0].begin() == 0 && subset[0].end() == 1 ) ;
				assert( parameters.size() == m_counts.size() ) ;
				m_parameters = parameters ;
			}
			
			double get_value_of_function() const {
				double result = 0.0 ;
				for( int i = 0; i < m_counts.size(); ++i ) {
					if( m_counts(i) > 0.0 ) {
						result += m_counts(i) * std::log( m_parameters( i )) ;
					}
				}
				
				return result ;
			}

			Vector get_value_of_first_derivative() const {
				Vector result = Vector::Zero( m_counts.size() ) ;
				for( int i = 0; i < result.size(); ++i ) {
					if( m_counts(i) > 0.0 ) {
						result(i) = m_counts(i) / m_parameters(i) ;
					}
				}
				return result ;
			}

			Matrix get_value_of_second_derivative() const {
				Matrix result = Matrix::Zero( m_counts.size(), m_counts.size() ) ;
				for( int i = 0; i < result.size(); ++i ) {
					if( m_counts(i) != 0.0 ) {
						result( i, i ) = -m_counts(i) / ( m_parameters(i) * m_parameters(i) ) ;
					}
				}
				return result ;
			}

			Vector get_MLE() const {
				return m_counts / m_counts.sum() ;
			}

			Vector parameters() const { return m_parameters ; }
			Vector const& get_data() const { return m_counts ; }
			
			std::string get_spec() const { return "Multinomial" ; }
		private:
			Vector m_counts ;
			Vector m_parameters ;
		} ;
	}
}

#endif
