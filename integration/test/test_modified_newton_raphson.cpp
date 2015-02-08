
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <boost/bind.hpp>
#include <boost/noncopyable.hpp>
#include "test_case.hpp"
#include <Eigen/Dense>
#include <Eigen/QR>
#include "integration/ModifiedNewtonRaphson.hpp"
#include "integration/AscentDirectionPicker.hpp"
#include "integration/ModifiedCholesky.hpp"


namespace impl {
	struct FunctionBase: public boost::noncopyable {
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
	} ;

	struct Multimodal1d_1: public FunctionBase {
		// This function f has two local maxima and one local minimum
		// f(x) = -x^4 + 2x^2
		// f'(x) = -4x^3 + 4x
		// f''(x) = -12 x^2 + 4
		// maxima are at \pm 1
		// minima is at 0
		void evaluate_at( Vector const& point, std::size_t number_of_derivatives = 0 ) {
			assert( point.size() == 1 ) ;
			m_point = point ;
		}
		double get_value_of_function() const {
			assert( m_point.size() == 1 ) ;
			return -std::pow( m_point(0), 4 ) + 2 * std::pow( m_point(0), 2 ) ;
		}
		Eigen::VectorXd get_value_of_first_derivative() const {
			assert( m_point.size() == 1 ) ;
			return Eigen::VectorXd::Constant( 1, -4.0 * std::pow( m_point(0), 3 ) + 4 * m_point(0) ) ; 
		}

		Eigen::MatrixXd get_value_of_second_derivative() const {
			return Eigen::MatrixXd::Constant(
				1, 1,
				-12.0 * std::pow( m_point(0), 2 ) + 4.0
			) ;
		}
	private:
		Vector m_point ;
	} ;

	struct Multimodal2d_1: public FunctionBase {
		// This function f has four local maxima and one local minimum.
		// We should walk to one of the maxima depending on where we start.
		// For appropriately modified NR, we should never walk to the minimum.
		// f' = -0.1 x^4 - 0.1 y^4 + x^2 + y^2
		// This has a minimum at (0,0)
		// And maxima at nonzero solutions of
		// -0.4 x^3 - 2x = 0 and -0.4 y^3 + 2y = 0
		// (x,y) = (\pm sqrt(5), \pm sqrt(5))
		//
		// First derivative is
		//
		// f' = ( -0.4 x^3 + 2x, -0.4 y^3 + 2y )
		//
		// Diagonals of second derivative are:
		// f'' = -1.2 x^2 + 2, -1.2 y^2 + 2
		//
		// With zeroes on off-diagonals
		//
		void evaluate_at( Vector const& point, std::size_t number_of_derivatives = 0 ) {
			assert( point.size() == 2 ) ;
			m_point = point ;
		}
		double get_value_of_function() const {
			return -0.1 * ( std::pow( m_point(0), 4 ) + std::pow( m_point(1), 4 ) )
				+ m_point(0) * m_point(0)
				+ m_point(1) * m_point(1)
			;
		}
		Eigen::VectorXd get_value_of_first_derivative() const {
			Eigen::VectorXd result ;
			result.setZero(2) ;
			result(0) = -0.4 * std::pow( m_point(0), 3 ) + 2.0 * m_point(0) ;
			result(1) = -0.4 * std::pow( m_point(1), 3 ) + 2.0 * m_point(1) ;
			return result ;
		}

		Eigen::MatrixXd get_value_of_second_derivative() const {
			Eigen::MatrixXd result ;
			result.setZero(2,2) ;
			result(0,0) = -1.2 * std::pow( m_point(0), 2 ) + 2.0 ;
			result(1,1) = -1.2 * std::pow( m_point(1), 2 ) + 2.0 ;
			return result ;
		}

	private:
		Vector m_point ;
	} ;
	
#if 0
	struct Rosenbrock: public FunctionBase {
		// Rosenbrock's function
		// has global minimum at (1,1)
		// f(x1,x2) = 100*(x2−x1^2)^2+(1−x1)^2
		// f'(x1,x2) = ( -4 * x1 * (x2−x1^2) -2 * x1,  )
		// f''(x1,x2) = ( -400*x1-2 
		void evaluate_at( Vector const& point, std::size_t number_of_derivatives = 0 ) {
			assert( point.size() == 2 ) ;
			m_point = point ;
		}
		double get_value_of_function() const {
			Vector const& x = m_point ;
			Vector x2 = m_point.array().square() ;
			return 100.0 * (x(1) - x2(0)) * (x(1) - x2(0)) + (1 - x(0)) * (1 - x(0) )  ;
		}

		Eigen::VectorXd get_value_of_first_derivative() const {
			assert(0) ;
		}

		Eigen::MatrixXd get_value_of_second_derivative() const {
			assert(0) ;
		}

	private:
		Vector m_point ;
	} ;
#endif
	
	template< typename Function >
	struct OvershootTarget: public boost::noncopyable {
		OvershootTarget( FunctionBase::Vector target ):
			m_target( target )
		{}
		
		FunctionBase::Vector compute( Function& function, FunctionBase::Vector const& v ) {
			// walk toward the target but go ten times too far.
			assert( v.size() == m_target.size() ) ;
			return 10.0 * ( m_target - v ) ;
		}
		
	private:
		FunctionBase::Vector m_target ;
	} ;
	
	template< typename Function >
	struct StoppingCondition {
		StoppingCondition( double tolerance  = 1E-12 ):
			m_tolerance( tolerance ),
			m_iterations(0)
			{}
		bool operator()(
			Function const& function,
			Eigen::VectorXd const& step
		) {
			++m_iterations ;
			// Stop when all elements of the derivative are close to 0
			// and the step size is also close to zero.
			return
				function.get_value_of_first_derivative().array().abs().maxCoeff() < m_tolerance
				&& step.array().abs().maxCoeff() < m_tolerance ;
		}
		std::size_t iterations() const { return m_iterations ; }
	private:
		double m_tolerance ;
		std::size_t m_iterations ;
	} ;
}

AUTO_TEST_CASE( test_newton_direction_1d ) {
	typedef impl::FunctionBase::Vector Vector ;
	typedef impl::FunctionBase::Matrix Matrix ;
	impl::Multimodal1d_1 function ;
	Eigen::ColPivHouseholderQR< Matrix > solver ;
	double const one_third = 0.3333333333333333333333333 ;
	Vector point( 1 ) ;
	for( point(0) = -10; point(0) < 10; point(0) += 0.01 ) {
		// Check that the solver always finds an ascent direction
		function.evaluate_at( point ) ;
		solver.compute( function.get_value_of_second_derivative() ) ;
		Vector h = solver.solve( -function.get_value_of_first_derivative() ) ;
		double directional_derivative = function.get_value_of_first_derivative().transpose() * h ;
		// f''(x) = -12 x^2 + 4
		// so 2nd derivative is negative for x < sqtr
		BOOST_CHECK(
			( std::abs( point(0) ) > std::sqrt( one_third ) && directional_derivative > 0 )
				|| 
			( std::abs( point(0) ) < std::sqrt( one_third ) && directional_derivative < 0 )
		) ;
	}
}

AUTO_TEST_CASE( test_ascent_direction_1d ) {
	typedef impl::FunctionBase::Vector Vector ;
	typedef impl::FunctionBase::Matrix Matrix ;
	impl::Multimodal1d_1 function ;
	integration::CholeskyOrEigenvalueSolver< impl::Multimodal1d_1 > solver( -0.1 ) ;

	Vector point( 1 ) ;
	for( point(0) = -10; point(0) < 10; point(0) += 0.01 ) {
		// Check that the solver always finds an ascent direction
		Vector h = solver.compute( function, point ) ;
		function.evaluate_at( point ) ;
		double directional_derivative = function.get_value_of_first_derivative().transpose() * h ;
		BOOST_CHECK_GT( directional_derivative, 0 ) ;
	}
}

namespace {
	template< typename Function >
	struct ModifiedCholeskySolver {
		typedef impl::FunctionBase::Vector Vector ;
		typedef impl::FunctionBase::Matrix Matrix ;

		Vector compute( Function& function, Vector const& point ) {
			function.evaluate_at( point ) ;
			m_solver.compute( -function.get_value_of_second_derivative() ) ;
			return m_solver.solve( function.get_value_of_first_derivative() ) ;
		} ;
		
	private:
		integration::ModifiedCholesky< Matrix > m_solver ;
	} ;
}

AUTO_TEST_CASE( test_modified_cholesky_direction_1d ) {
	typedef impl::FunctionBase::Vector Vector ;
	typedef impl::FunctionBase::Matrix Matrix ;
	impl::Multimodal1d_1 function ;
	ModifiedCholeskySolver< impl::Multimodal1d_1 > solver ;

	Vector point( 1 ) ;
	for( point(0) = -10; point(0) < 10; point(0) += 0.01 ) {
		// Check that the solver always finds an ascent direction
		Vector h = solver.compute( function, point ) ;
		double directional_derivative = function.get_value_of_first_derivative().transpose() * h ;
		BOOST_CHECK_GT( directional_derivative, 0 ) ;
	}
}

AUTO_TEST_CASE( test_ascent_direction_2d ) {
	typedef impl::FunctionBase::Vector Vector ;
	typedef impl::FunctionBase::Matrix Matrix ;
	impl::Multimodal2d_1 function ;
	integration::CholeskyOrEigenvalueSolver< impl::Multimodal2d_1 > solver( -0.1 ) ;

	Vector point( 2 ) ;
	for( point(0) = -10; point(0) < 10; point(0) += 0.025 ) {
		for( point(1) = -10; point(1) < 10; point(1) += 0.025 ) {
			Vector h = solver.compute( function, point ) ;
			function.evaluate_at( point ) ;
			double const directional_derivative = function.get_value_of_first_derivative().transpose() * h ;
			BOOST_CHECK_GT( directional_derivative, 0 ) ;
		}
	}
}

AUTO_TEST_CASE( test_modified_cholesky_direction_2d ) {
	typedef impl::FunctionBase::Vector Vector ;
	typedef impl::FunctionBase::Matrix Matrix ;
	impl::Multimodal2d_1 function ;
	ModifiedCholeskySolver< impl::Multimodal2d_1 > solver ;

	Vector point( 2 ) ;
	for( point(0) = -10; point(0) < 10; point(0) += 0.025 ) {
		for( point(1) = -10; point(1) < 10; point(1) += 0.025 ) {
			Vector h = solver.compute( function, point ) ;
			function.evaluate_at( point ) ;
			double const directional_derivative = function.get_value_of_first_derivative().transpose() * h ;
			BOOST_CHECK_GT( directional_derivative, 0 ) ;
		}
	}
}

AUTO_TEST_CASE( test_ascent_direction_newton_raphson_overshoot_1d ) {
	typedef impl::FunctionBase::Vector Vector ;
	impl::Multimodal1d_1 function ;
	impl::StoppingCondition< impl::Multimodal1d_1 > stoppingCondition( 1E-12 ) ;
	integration::CholeskyOrEigenvalueSolver< impl::Multimodal1d_1 > solver ;

	// For 1d function 1, maxima are at -1 and 1.
	Vector target ;
	target.setZero(1) ;
	target(0) = -1 ;

	{
		impl::OvershootTarget< impl::Multimodal1d_1 > overshoot( target ) ;
		Vector start = Vector::Constant( 1, -0.001 ) ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			overshoot,
			stoppingCondition
		) ;

		BOOST_CHECK_CLOSE( result(0), -1, 1E-7 ) ;
	}
}

AUTO_TEST_CASE( test_modified_cholesky_newton_raphson_overshoot_1d ) {
	typedef impl::FunctionBase::Vector Vector ;
	impl::Multimodal1d_1 function ;
	impl::StoppingCondition< impl::Multimodal1d_1 > stoppingCondition( 1E-12 ) ;
	ModifiedCholeskySolver< impl::Multimodal1d_1 > solver ;

	// For 1d function 1, maxima are at -1 and 1.
	Vector target ;
	target.setZero(1) ;
	target(0) = -1 ;

	{
		impl::OvershootTarget< impl::Multimodal1d_1 > overshoot( target ) ;
		Vector start = Vector::Constant( 1, -0.001 ) ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			overshoot,
			stoppingCondition
		) ;

		BOOST_CHECK_CLOSE( result(0), -1, 1E-7 ) ;
	}
}

AUTO_TEST_CASE( test_ascent_direction_newton_raphson_overshoot_2d ) {
	typedef impl::FunctionBase::Vector Vector ;
	impl::Multimodal2d_1 function ;
	impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
	integration::CholeskyOrEigenvalueSolver< impl::Multimodal2d_1 > solver ;

	Vector maximum ;
	maximum.setZero(2) ;
	maximum(0) = -sqrt(5) ;
	maximum(1) = -sqrt(5) ;

	{
		impl::OvershootTarget< impl::Multimodal2d_1 > overshoot( maximum ) ;
		Vector start = Vector::Constant(2,-1) ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			overshoot,
			stoppingCondition
		) ;

		BOOST_CHECK_CLOSE( result(0), -std::sqrt(5), 1E-5 ) ;
		BOOST_CHECK_CLOSE( result(1), -std::sqrt(5), 1E-5 ) ;
	}
}

AUTO_TEST_CASE( test_modified_cholesky_newton_raphson_overshoot_2d ) {
	typedef impl::FunctionBase::Vector Vector ;
	impl::Multimodal2d_1 function ;
	impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
	ModifiedCholeskySolver< impl::Multimodal2d_1 > solver ;

	Vector maximum ;
	maximum.setZero(2) ;
	maximum(0) = -sqrt(5) ;
	maximum(1) = -sqrt(5) ;

	{
		impl::OvershootTarget< impl::Multimodal2d_1 > overshoot( maximum ) ;
		Vector start = Vector::Constant(2,-1) ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			overshoot,
			stoppingCondition
		) ;

		BOOST_CHECK_CLOSE( result(0), -std::sqrt(5), 1E-5 ) ;
		BOOST_CHECK_CLOSE( result(1), -std::sqrt(5), 1E-5 ) ;
	}
}

AUTO_TEST_CASE( test_ascent_direction_newton_raphson_2d ) {
	typedef impl::FunctionBase::Vector Vector ;
	impl::Multimodal2d_1 function ;
	impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
	integration::CholeskyOrEigenvalueSolver< impl::Multimodal2d_1 > solver ;

	{		
		impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
		Vector start = Vector::Constant(2,-1) ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			solver,
			stoppingCondition
		) ;
			
		BOOST_CHECK_CLOSE( result(0), -std::sqrt(5), 1E-5 ) ;
		BOOST_CHECK_CLOSE( result(1), -std::sqrt(5), 1E-5 ) ;
		std::cerr << "test_modified_cholesky_newton_raphson_2d: took " << stoppingCondition.iterations() << " iterations.\n" ;
	}
	{	
		impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
		Vector start = Vector::Constant(2,-1) ;
		start(1) = -0.1 ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			solver,
			stoppingCondition
		) ;
			
		BOOST_CHECK_CLOSE( result(0), -std::sqrt(5), 1E-5 ) ;
		BOOST_CHECK_CLOSE( result(1), -std::sqrt(5), 1E-5 ) ;
		std::cerr << "test_modified_cholesky_newton_raphson_2d: took " << stoppingCondition.iterations() << " iterations.\n" ;
	}
	{
		impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
		Vector start = Vector::Constant(2,-1) ;
		start(0) = -100 ;
		start(1) = -0.1 ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			solver,
			stoppingCondition
		) ;
			
		BOOST_CHECK_CLOSE( result(0), -std::sqrt(5), 1E-5 ) ;
		BOOST_CHECK_CLOSE( result(1), -std::sqrt(5), 1E-5 ) ;
		std::cerr << "test_modified_cholesky_newton_raphson_2d: took " << stoppingCondition.iterations() << " iterations.\n" ;
	}
}

AUTO_TEST_CASE( test_modified_cholesky_newton_raphson_2d ) {
	typedef impl::FunctionBase::Vector Vector ;
	impl::Multimodal2d_1 function ;
	ModifiedCholeskySolver< impl::Multimodal2d_1 > solver ;

	{		
		impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
		Vector start = Vector::Constant(2,-1) ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			solver,
			stoppingCondition
		) ;
			
		BOOST_CHECK_CLOSE( result(0), -std::sqrt(5), 1E-5 ) ;
		BOOST_CHECK_CLOSE( result(1), -std::sqrt(5), 1E-5 ) ;
		std::cerr << "test_modified_cholesky_newton_raphson_2d: took " << stoppingCondition.iterations() << " iterations.\n" ;
	}
	{	
		impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
		Vector start = Vector::Constant(2,-1) ;
		start(1) = -0.1 ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			solver,
			stoppingCondition
		) ;
			
		BOOST_CHECK_CLOSE( result(0), -std::sqrt(5), 1E-5 ) ;
		BOOST_CHECK_CLOSE( result(1), -std::sqrt(5), 1E-5 ) ;
		std::cerr << "test_modified_cholesky_newton_raphson_2d: took " << stoppingCondition.iterations() << " iterations.\n" ;
	}
	{
		impl::StoppingCondition< impl::Multimodal2d_1 >  stoppingCondition( 1E-12 ) ;
		Vector start = Vector::Constant(2,-1) ;
		start(0) = -100 ;
		start(1) = -0.1 ;
		Vector result = integration::find_maximum_by_modified_newton_raphson_with_line_search(
			function,
			start,
			solver,
			stoppingCondition
		) ;
			
		BOOST_CHECK_CLOSE( result(0), -std::sqrt(5), 1E-5 ) ;
		BOOST_CHECK_CLOSE( result(1), -std::sqrt(5), 1E-5 ) ;
		std::cerr << "test_modified_cholesky_newton_raphson_2d: took " << stoppingCondition.iterations() << " iterations.\n" ;
	}
}

ModifiedCholeskySolver< impl::Multimodal2d_1 > solver ;