#include <iostream>
#include <boost/bind.hpp>
#include "test_case.hpp"
#include "integration/NewtonRaphson.hpp"
#include <Eigen/Dense>

namespace impl {
	double function_1( double v ) {
		return -v*v + v + 2.0 ;
	}

	double derivative_1( double v ) {
		return -2.0 * v + 1 ;
	}
	
	double quadratic( double v, double a, double b, double c ) {
		return a*v*v + b*v + c ;
	}
	
	double quadratic_derivative( double v, double a, double b ) {
		return 2*a*v + b ;
	}
	
	Eigen::VectorXd affine_map( Eigen::VectorXd const& point, Eigen::MatrixXd const& matrix, Eigen::VectorXd const& offset ) {
		return ( matrix * ( point - offset ) ) ;
	}

	Eigen::MatrixXd affine_derivative( Eigen::VectorXd const& point, Eigen::MatrixXd const& matrix ) {
		return matrix ;
	}

	Eigen::VectorXd another_map( Eigen::VectorXd const& point, Eigen::MatrixXd const& matrix, Eigen::VectorXd const& offset ) {
		// A linear map plus quadratic deviation.
		Eigen::VectorXd v( 2 ) ;
		v << 	-( ( point - offset )( 0 ) * ( point - offset )( 0 ) ),
				-( ( point - offset )( 1 ) * ( point - offset )( 1 ) ) ;
		return ( matrix * ( point - offset ) ) + v ;
	}

	Eigen::MatrixXd another_derivative( Eigen::VectorXd const& point, Eigen::MatrixXd const& matrix, Eigen::VectorXd const& offset ) {
		Eigen::MatrixXd M( 2, 2 ) ;
		M <<	-2 * ( point( 0 ) - offset( 0 ) ),	0,
				0,									-2 * ( point( 1 ) - offset( 1 ) ) ;
		return matrix + M ;
	}
}

void test_newton_raphson_1d() {
	std::cerr << "test_newton_raphson_1d(): finding roots of  5 x - 2...\n" ;
	for( double epsilon = 0.1; epsilon > 0.0000000000001; epsilon /= 10 ) {
		std::cerr << "epsilon = " << epsilon << ".\n" ;
		double root = integration::find_root_by_newton_raphson(
			boost::bind( impl::quadratic, _1, 0, 5, -2 ),
			boost::bind( impl::quadratic_derivative, _1, 0, 5 ),
			-10.0,
			epsilon
		) ;
		
		std::cerr << "root: " << root << ".\n" ;
		TEST_ASSERT( std::abs( root - ( 2.0 / 5.0 )) < epsilon ) ;
	}
	
	std::cerr << "test_newton_raphson_1d(): finding roots of 2 + x - x^2...\n" ;
	for( double epsilon = 0.1; epsilon > 0.0000000000001; epsilon /= 10 ) {
		std::cerr << "epsilon = " << epsilon << ".\n" ;
		double left_root = integration::find_root_by_newton_raphson(
			boost::bind( impl::quadratic, _1, -1, 1, 2 ),
			boost::bind( impl::quadratic_derivative, _1, -1, 1 ),
			-10.0,
			epsilon
		) ;
		std::cerr << "Left root: " << left_root << ".\n" ;
		TEST_ASSERT( std::abs( left_root + 1 ) < epsilon ) ;

		double right_root = integration::find_root_by_newton_raphson(
			&impl::function_1,
			&impl::derivative_1,
			10.0,
			epsilon
		) ;

		std::cerr << "Right root: " << right_root << ".\n" ;
		TEST_ASSERT( std::abs( right_root - 2 ) < epsilon ) ;
	}
}

void test_newton_raphson_nd_1() {
	std::cerr << "test_newton_raphson_nd_1(): finding roots of \\Sigma x = 0...\n" ;
	
	Eigen::MatrixXd sigma ;
	sigma.resize( 2, 2 ) ;
	sigma( 0, 0 ) = 5 ;
	sigma( 1, 1 ) = 5 ;
	sigma( 0, 1 ) = 0.5 ;
	sigma( 1, 0 ) = 0.5 ;

	Eigen::VectorXd actual_root = Eigen::VectorXd( 2 ) ;
	actual_root << 50, 50 ;
	
	Eigen::VectorXd initial_point( 2 ) ;
	initial_point << -100.0, -100.0 ;
	for( double epsilon = 0.1; epsilon > 0.0000000000001; epsilon /= 10 ) {
		std::cerr << "==================== epsilon = " << epsilon << ".\n" ;
		Eigen::VectorXd root = integration::find_root_by_newton_raphson(
			boost::bind( impl::affine_map, _1, sigma, actual_root ),
			boost::bind( impl::affine_derivative, _1, sigma ),
			initial_point,
			epsilon
		) ;
	
		std::cerr << "root is " << root << ".\n" ;
		TEST_ASSERT( ( root - actual_root ).norm() < epsilon ) ;
	}
}

void test_newton_raphson_nd_2() {
	std::cerr << "test_newton_raphson_nd_2(): finding roots of \\Sigma x = 0...\n" ;
	
	Eigen::MatrixXd sigma ; // variance-covariance matrix
	sigma.resize( 2, 2 ) ;
	sigma( 0, 0 ) = 5 ;
	sigma( 1, 1 ) = 5 ;
	sigma( 0, 1 ) = 0 ;
	sigma( 1, 0 ) = 0 ;

	Eigen::VectorXd actual_root = Eigen::VectorXd( 2 ) ;
	actual_root << 50, 50 ;
	
	Eigen::VectorXd initial_point( 2 ) ;
	initial_point << -100.0, -100.0 ;
	for( double epsilon = 0.1; epsilon > 0.000000000000001; epsilon /= 10 ) {
		std::cerr << "============== epsilon = " << epsilon << ".\n" ;
		Eigen::VectorXd root = integration::find_root_by_newton_raphson(
			boost::bind( impl::another_map, _1, sigma, actual_root ),
			boost::bind( impl::another_derivative, _1, sigma, actual_root ),
			initial_point,
			epsilon
		) ;
	
		std::cerr << "root is " << root << ".\n" ;
		TEST_ASSERT( impl::another_map( root, sigma, actual_root ).maxCoeff() < epsilon ) ;
	}
}

int main(){
	test_newton_raphson_1d() ;
	test_newton_raphson_nd_1() ;
	test_newton_raphson_nd_2() ;
}