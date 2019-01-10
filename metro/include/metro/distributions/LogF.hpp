
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_DISTRIBUTIONS_LOGF_HPP
#define METRO_DISTRIBUTIONS_LOGF_HPP

#include "metro/UnivariateLogDensity.hpp"

namespace metro {
	namespace distributions {
		struct LogF: public UnivariateLogDensity {
		public:
			enum Normalisation { ePDF, eZeroAtMean } ;
			static UniquePtr create( double nu1, double nu2 ) ;
		public:
			LogF( double nu1, double nu2 ) ;
			LogF( LogF const& other ) ;
			LogF& operator=( LogF const& other ) ;

		public:
			void evaluate_at( double x ) ;

			double get_value_of_function() const {
				return m_value_of_function ;
			}

			double get_value_of_first_derivative() const {
				return m_value_of_first_derivative ;
			}

			double get_value_of_second_derivative() const {
				return m_value_of_second_derivative ;
			}
			
			std::string get_summary() const ;
				
		private:
			double m_alpha ;
			double m_beta ;
			double m_constant ;
			double m_value_of_function ;
			double m_value_of_first_derivative ;
			double m_value_of_second_derivative ;
		} ;
	}
}

#endif
