
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_STOPPINGCONDITION_HPP
#define METRO_STOPPINGCONDITION_HPP

namespace metro {
	// Stop when repeat values differ by no more than the given amount
	// Or when a maximum number of calls is made.
	struct ValueStabilisesStoppingCondition {
		ValueStabilisesStoppingCondition( double const tolerance, std::size_t max_iterations = 10000 ):
			m_tolerance( tolerance ),
			m_max_iterations( max_iterations )
		{
			reset() ;
		}
			
		// Return
		bool operator()( double value ) {
			m_converged = std::abs( value - m_value ) < m_tolerance ;
			m_value = value ;

			return(
				m_converged || ( ++m_iteration >= m_max_iterations ) || ( value != value )
			) ;
		}
	
		bool converged() const {
			return m_converged ;
		}
	
		std::size_t iterations() const {
			return m_iteration ;
		}

		void reset() {
			m_converged = false ;
			m_iteration = 0 ;
			m_value = std::numeric_limits< double >::infinity() ;
		}

	private:
		double const m_tolerance ;
		std::size_t const m_max_iterations ;

		std::size_t m_iteration ;
		double m_value ;
		bool m_converged ;
	} ;
}

#endif

