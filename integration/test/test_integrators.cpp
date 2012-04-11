
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include "test_case.hpp"
#include "fputils/floating_point_utils.hpp"
#include "integration/Integrator.hpp"
#include "integration/GanderGautschiAdaptiveIntegrator.hpp"
#include "integration/GanderGautschiAdaptiveLogIntegrator.hpp"


namespace integrands {
	// Some functions to integrate.
	double exponential( double a ) { return std::exp( a ) ; }

	double log_of_exponential( double a ) { return a ; }

	struct Exponential
	{
		Exponential( double scale ): m_scale( scale ), m_number_of_evaluations(0) {}
		double operator()( double a ) const { ++m_number_of_evaluations; return std::exp( m_scale * a ) ; }
		~Exponential() { std::cerr << "Exponential( " << m_scale << "): evaluated " << m_number_of_evaluations << " times.\n" ; }
	private:
		double const m_scale ;
		mutable std::size_t m_number_of_evaluations ;
	} ;

	struct LogOfExponential
	{
		LogOfExponential( double scale ): m_scale( scale ), m_number_of_evaluations(0) {}
		double operator()( double a ) const { ++m_number_of_evaluations ; return m_scale * a ; }
		~LogOfExponential() { std::cerr << "LogOfExponential( " << m_scale << "): evaluated " << m_number_of_evaluations << " times.\n" ; }
	private:
		double const m_scale ;
		mutable std::size_t m_number_of_evaluations ;
	} ;
	
	double square_root( double a ) { return std::sqrt(a) ; }
	double log_of_square_root( double a ) { return std::log( std::sqrt(a )) ; }
}

AUTO_TEST_CASE( test_gander_gauschi )
{
	std::cerr << "Test Gander-Gauschi integrator...\n" ;

	std::vector< double > tolerances ;
	tolerances.push_back( 1E-4 ) ;
	tolerances.push_back( 1E-6 ) ;
	tolerances.push_back( 1E-8 ) ;
	tolerances.push_back( 1E-10 ) ;
	tolerances.push_back( 10.0 * std::numeric_limits< double >::epsilon() ) ;

	for( std::size_t i = 0; i < tolerances.size(); ++i ) {
		double const& tolerance = tolerances[i] ;
		std::cerr << "test_gander_gauschi: using tolerance " << std::fixed << std::setprecision(20) << tolerance << "...\n" ;
		integration::GanderGautschiAdaptiveIntegrator integrator( tolerance ) ;

		// We do not always get within tolerance.  But we do always get within 10 times it.
		// This corresponds to grey bars in the histograms on Gander - Gautschi's paper.
		// Note that MATLAB has a similar implementation with similar results.
		double global_error = tolerance * 10 ;

		// Test exponential integration
		{
			double result = integrator( &integrands::exponential, 0.0, 1.0 ) ;
			double expected = std::exp(1) - 1.0 ;
			if( std::abs( expected - result ) > tolerance ) {
				std::cerr << "Warning: missed tolerance for \\int_0^1 e^{x} by " << std::abs( expected - result ) << ": expected " << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < global_error ) ;
			}
		}
		
		{
			double result = integrator( integrands::Exponential( 5.0 ), 0.0, 1.0 ) ;
			double expected = ( std::exp( 5.0 ) - 1.0 ) / 5.0 ;
			if( std::abs( expected - result ) > tolerance ) {
				std::cerr << "Warning: missed tolerance for \\int_0^1 e^{5x} by " << std::abs( expected - result ) << ": expected " << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < global_error ) ;
			}
		}

		try {
			integrator( integrands::Exponential( 1000.0 ), 0.0, 1.0 ) ;
			TEST_ASSERT( 0 ) ;
		}
		catch( integration::IntervalTooSmallToSubdivideError const& e ) {
			// ok, this should fail.
		}

		{
			double result = integrator( &integrands::square_root, 0.0, 1.0 ) ;
			double expected = 2.0 / 3.0 ;
			if( std::abs( expected - result ) > tolerance ) {			
				std::cerr << "Warning: missed tolerance for \\int_0^1 \\sqrt(x) by " << std::abs( expected - result ) << ": expected " << std::fixed << std::setprecision(20) << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < global_error ) ;
			}
		}

		{
			double result = integrator( &integrands::square_root, 0.0, 10.0 ) ;
			double expected = 2.0 * std::pow( 10, 3.0 / 2.0 ) / 3.0 ;
			if( std::abs( expected - result ) > tolerance ) {			
				std::cerr << "Warning: missed tolerance for \\int_0^10 \\sqrt(x) by " << std::abs( expected - result ) << ": expected " << std::fixed << std::setprecision(20) << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < (10 * global_error) ) ;
			}
		}
	}
}


AUTO_TEST_CASE( test_gander_gauschi_log )
{
	std::cerr << "Test Gander-Gauschi log-integrator...\n" ;

	std::vector< double > tolerances ;
	tolerances.push_back( 1E-4 ) ;
	tolerances.push_back( 1E-6 ) ;
	tolerances.push_back( 1E-8 ) ;
	tolerances.push_back( 1E-10 ) ;
	tolerances.push_back( 10.0 * std::numeric_limits< double >::epsilon() ) ;

	for( std::size_t i = 0; i < tolerances.size(); ++i ) {
		double const& tolerance = tolerances[i] ;
		std::cerr << "test_gander_gauschi: using tolerance " << std::fixed << std::setprecision(20) << tolerance << "...\n" ;
		integration::GanderGautschiAdaptiveLogIntegrator integrator( tolerance ) ;

		// We do not always get within tolerance.  But we do always get within 10 times it.
		// This corresponds to grey bars in the histograms on Gander - Gautschi's paper.
		// Note that MATLAB has a similar implementation with similar results.
		double global_error = tolerance * 10 ;

		// Test exponential integration
		{
			double result = integrator( &integrands::log_of_exponential, 0.0, 1.0 ) ;
			double expected = std::log( std::exp(1) - 1.0 ) ;
			if( std::abs( expected - result ) > tolerance ) {
				std::cerr << "Warning: missed tolerance for \\int_0^1 e^{x} by " << std::abs( expected - result ) << ": expected " << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < global_error ) ;
			}
		}
		
		{
			double result = integrator( integrands::LogOfExponential( 5.0 ), 0.0, 1.0 ) ;
			double expected = std::log( ( std::exp( 5.0 ) - 1.0 ) ) - std::log( 5.0 ) ;
			if( std::abs( expected - result ) > tolerance ) {
				std::cerr << "Warning: missed tolerance for \\int_0^1 e^{5x} by " << std::abs( expected - result ) << ": expected " << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < global_error ) ;
			}
		}
		
		{
			double result = integrator( integrands::LogOfExponential( 1000.0 ), 0.0, 1.0 ) ;
			double expected = fputils::log_diff_exp( 1000.0, 0.0 ) - std::log( 1000.0 ) ;
			if( std::abs( expected - result ) > tolerance ) {
				std::cerr << "Warning: missed tolerance for \\int_0^1 e^{1000x} by " << std::abs( expected - result ) << ": expected " << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < global_error ) ;
			}
		}

		{
			double result = integrator( &integrands::log_of_square_root, 0.0, 1.0 ) ;
			double expected = std::log( 2.0 ) - std::log( 3.0 ) ;
			if( std::abs( expected - result ) > tolerance ) {			
				std::cerr << "Warning: missed tolerance for \\int_0^1 \\sqrt(x) by " << std::abs( expected - result ) << ": expected " << std::fixed << std::setprecision(20) << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < global_error ) ;
			}
		}

		{
			double result = integrator( &integrands::log_of_square_root, 0.0, 10.0 ) ;
			double expected = std::log( 2.0 ) + std::log( std::pow( 10, 3.0 / 2.0 ) ) - std::log( 3.0 ) ;
			if( std::abs( expected - result ) > tolerance ) {			
				std::cerr << "Warning: missed tolerance for \\int_0^10 \\sqrt(x) by " << std::abs( expected - result ) << ": expected " << std::fixed << std::setprecision(20) << expected << ", got: " << result << ".\n" ;
				TEST_ASSERT( std::abs( expected - result ) < (10 * global_error) ) ;
			}
		}
	}
}

