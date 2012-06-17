
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_SUM_HPP
#define METRO_LIKELIHOOD_SUM_HPP

#include <memory>
#include <Eigen/Core>
#include <boost/ptr_container/ptr_vector.hpp>

namespace metro {
	namespace likelihood {
		template< typename Scalar, typename Vector, typename Matrix >
		struct Sum: public metro::LogLikelihood< Scalar, Vector, Matrix >
		{
			typedef metro::LogLikelihood< Scalar, Vector, Matrix > LogLikelihood ;
			void add_term( LogLikelihood::UniquePtr term ) ;
			
			Scalar get_value_of_function() const {
				Scalar result = 0.0 ;
				for( std::size_t i = 0; i < m_terms.size(); ++i ) {
					result += m_terms[i].get_value_of_function() ;
				}
				return result ;
			}

			Vector get_value_of_first_derivative() const {
				Vector result ;
				for( std::size_t i = 0; i < m_terms.size(); ++i ) {
					if( i == 0 ) {
						result = m_terms[i].get_value_of_first_derivative() ;
					} else {
						result += m_terms[i].get_value_of_first_derivative() ;
					}
				}
				return result ;
			}

			Matrix get_value_of_second_derivative() const {
				Matrix result ;
				for( std::size_t i = 0; i < m_terms.size(); ++i ) {
					if( i == 0 ) {
						result = m_terms[i].get_value_of_second_derivative() ;
					} else {
						result += m_terms[i].get_value_of_second_derivative() ;
					}
				}
				return result ;
			}
			
			std::string get_spec() const {
				std::string result =  "Sum( " ;
				for( std::size_t i = 0; i < m_terms.size(); ++i ) {
					if( i > 0 ) {
						result += ", " ;
					}
					result += m_terms[i].get_spec() ;
				}
				return result + ")" ;
			}

		private:
			boost::ptr_vector< LogLikelihood > m_terms ;
		} ;
	}
}

#endif
