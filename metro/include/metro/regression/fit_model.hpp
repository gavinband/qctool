
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_REGRESSION_FIT_MODEL_HPP
#define METRO_REGRESSION_FIT_MODEL_HPP

#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "boost/function.hpp"
#include "metro/ModifiedCholesky.hpp"
#include "metro/SmoothFunction.hpp"

namespace metro {
	namespace regression {
		struct Stepper {
		public:
			typedef SmoothFunction Function ;
			typedef Eigen::VectorXd Vector ;
			typedef Eigen::MatrixXd Matrix ;
			typedef boost::function< void( int iteration, double ll, double target_ll, Vector const& point, Vector const& derivative, Vector const& step, bool converged ) > Tracer ;
		public:
			virtual ~Stepper() ;
			virtual bool diverged() const = 0 ;
			virtual bool step( Function& function, Vector const& point, Vector* result ) = 0 ;
			virtual std::size_t number_of_iterations() const = 0 ;
		} ;

		struct CholeskyStepper: public Stepper {
		public:
			CholeskyStepper( double tolerance, int max_iterations, Tracer tracer ) ;

			bool diverged() const ;
			std::size_t number_of_iterations() const ;

			bool step( Function& function, Vector const& point, Vector* result ) ;
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

		std::pair< bool, int > fit_model(
			metro::SmoothFunction& ll,
			std::string const& model_name,
			Eigen::VectorXd const& starting_point,
			Stepper& stopping_condition,
			std::vector< std::string >* comments
		) ;
	}
}

#endif
