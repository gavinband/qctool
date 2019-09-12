
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_CHOLESKY_STEPPER_HPP
#define METRO_CHOLESKY_STEPPER_HPP

#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "boost/function.hpp"
#include "metro/ModifiedCholesky.hpp"
#include "metro/SmoothFunction.hpp"
#include "metro/fit_model.hpp"

namespace metro {
	struct CholeskyStepper: public Stepper {
	public:
		CholeskyStepper( double tolerance, int max_iterations, Tracer tracer = Tracer()) ;

		bool step( Function& function, Vector const& point, Vector* result ) ;

		bool diverged() const ;
		std::size_t number_of_iterations() const ;
		void reset() ;

	private:
		double const m_tolerance ;
		int const m_max_iterations ;
		Tracer m_tracer ;
		int m_iteration ;
		metro::ModifiedCholesky< Matrix > m_solver ;
		//Eigen::LDLT< Matrix > m_solver ;
		//Eigen::ColPivHouseholderQR< Matrix > m_solver ;
		double m_target_ll ;
	} ;
}

#endif
