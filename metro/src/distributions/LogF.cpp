
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <cmath>
#include <cassert>
#include <boost/format.hpp>
#include "metro/distributions/LogF.hpp"

namespace metro {
	namespace distributions {
		LogF::UniquePtr LogF::create( double nu1, double nu2 ) {
			return LogF::UniquePtr( new LogF( nu1, nu2 )) ;
		}

		LogF::LogF( double nu1, double nu2 ):
			m_alpha( nu1 / 2 ),
			m_beta( nu2 / 2 ),
			m_constant(0)
		{
			m_constant = -(m_alpha + m_beta) * std::log(0.5) ;
		}
		
		LogF::LogF( LogF const& other ):
			m_alpha( other.m_alpha ),
			m_beta( other.m_beta ),
			m_constant( other.m_constant )
		{}
			
		LogF& LogF::operator=( LogF const& other ) {
			m_alpha = other.m_alpha ;
			m_beta = other.m_beta ;
			m_constant = other.m_constant ;
			return *this ;
		}
		

		void LogF::evaluate_at( double x ) {
			double const l = 1.0 / (1.0 + std::exp(-x)) ;
			m_value_of_function = m_constant + ( m_alpha * std::log(l) + m_beta * std::log( 1-l ) );
			m_value_of_first_derivative = m_alpha - (m_alpha + m_beta) * l ;
			m_value_of_second_derivative = -(m_alpha + m_beta) * l * (1-l) ;
		}
		
		std::string LogF::get_summary() const {
			return(
				boost::format( "logF( %.5d, %.5d )" ) % (m_alpha*2) % (m_beta*2)
			).str() ;
		}
	}
}
