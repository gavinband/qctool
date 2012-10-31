
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <Eigen/Core>
#include "test_case.hpp"
#include "metro/mean_and_covariance.hpp"

AUTO_TEST_CASE( test_mean_and_covariance_with_missingness ) {
	using metro::compute_mean_and_covariance ;

	double const tolerance = 0.000000000001 ;
	Eigen::RowVectorXd mean ;
	Eigen::MatrixXd covariance ;
	Eigen::MatrixXd X ;
	Eigen::MatrixXd nonmissingness ;

	{
		X.resize( 1, 2 ) ;
		X << 1, 2 ;
		nonmissingness.resize( X.rows(), X.cols() ) ;
		nonmissingness.setOnes() ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( mean(0), 1 ) ;
		BOOST_CHECK_EQUAL( mean(1), 2 ) ;
		BOOST_CHECK( std::isnan( covariance(0,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,1) ) ) ;
	
		nonmissingness(0) = 0 ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK( std::isnan( mean(0) ) ) ;
		BOOST_CHECK_EQUAL( mean(1), 2 ) ;
		BOOST_CHECK( std::isnan( covariance(0,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,1) ) ) ;
	}
	
	{
		X.resize( 2, 2 ) ;
		X <<	0, 1,
				1, 2 ;
		nonmissingness.resize( X.rows(), X.cols() ) ;
		nonmissingness.setOnes() ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean(0), 0.5 ) ;
		BOOST_CHECK_EQUAL( mean(1), 1.5 ) ;
		BOOST_CHECK_EQUAL( covariance(0,0), 0.5 ) ;
		BOOST_CHECK_EQUAL( covariance(0,1), 0.5 ) ;
		BOOST_CHECK_EQUAL( covariance(1,1), 0.5 ) ;
		BOOST_CHECK_EQUAL( covariance(0,1), covariance(1,0) ) ;

		nonmissingness(1,0) = 0 ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean(0), 0 ) ;
		BOOST_CHECK_EQUAL( mean(1), 1.5 ) ;
		BOOST_CHECK( std::isnan( covariance(0,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
		BOOST_CHECK_EQUAL( covariance(1,1), 0.5 ) ;

		nonmissingness(0,1) = 0 ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean(0), 0 ) ;
		BOOST_CHECK_EQUAL( mean(1), 2 ) ;
		BOOST_CHECK( std::isnan( covariance(0,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,1) ) ) ;

		nonmissingness(0,0) = 0 ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK( std::isnan( mean(0) ) ) ;
		BOOST_CHECK_EQUAL( mean(1), 2 ) ;
		BOOST_CHECK( std::isnan( covariance(0,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,1) ) ) ;

		nonmissingness(1,1) = 0 ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK( std::isnan( mean(0) ) ) ;
		BOOST_CHECK( std::isnan( mean(1) ) ) ;
		BOOST_CHECK( std::isnan( covariance(0,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,1) ) ) ;
	}
	
	{
		X.resize( 2, 2 ) ;
		X <<	0,		1,
				0.5,	2 ;
		nonmissingness.resize( X.rows(), X.cols() ) ;
		nonmissingness.setOnes() ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean(0), 0.25 ) ;
		BOOST_CHECK_EQUAL( mean(1), 1.5 ) ;
		BOOST_CHECK_EQUAL( covariance(0,0), 0.125 ) ;
		BOOST_CHECK_EQUAL( covariance(0,1), 0.25 ) ;
		BOOST_CHECK_EQUAL( covariance(1,1), 0.5 ) ;
		BOOST_CHECK_EQUAL( covariance(0,1), covariance(1,0) ) ;

		compute_mean_and_covariance( X.transpose(), nonmissingness.transpose(), mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean(0), 0.5 ) ;
		BOOST_CHECK_EQUAL( mean(1), 1.25 ) ;
		BOOST_CHECK_EQUAL( covariance(0,0), 0.5 ) ;
		BOOST_CHECK_EQUAL( covariance(0,1), 0.75 ) ;
		BOOST_CHECK_EQUAL( covariance(1,1), 1.125 ) ;
		BOOST_CHECK_EQUAL( covariance(0,1), covariance(1,0) ) ;
		
		nonmissingness(1,0) = 0 ;
		compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean(0), 0 ) ;
		BOOST_CHECK_EQUAL( mean(1), 1.5 ) ;
		BOOST_CHECK( std::isnan( covariance(0,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
		BOOST_CHECK_EQUAL( covariance(1,1), 0.5 ) ;

		compute_mean_and_covariance( X.transpose(), nonmissingness.transpose(), mean, covariance ) ;
		BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
		BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
		BOOST_CHECK_EQUAL( mean(0), 0.5 ) ;
		BOOST_CHECK_EQUAL( mean(1), 2 ) ;
		BOOST_CHECK_EQUAL( covariance(0,0), 0.5 ) ;
		BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
		BOOST_CHECK( std::isnan( covariance(1,1) ) ) ;
	}
	
	{
		X.resize( 20, 2 ) ;
		X <<
			-0.3105243,9.721345,
			0.9230904,10.76971,
			-0.2534216,9.653925,
			-1.579065,8.344736,
			1.817072,11.8996,
			0.4140253,10.27366,
			0.606559,10.52624,
			1.750472,11.82706,
			-2.246558,7.49725,
			0.8391379,10.90392,
			-0.5266163,9.140894,
			1.020504,11.02627,
			-0.6118804,9.614296,
			0.5836782,10.69602,
			-1.59382,8.367049,
			-0.1051234,9.739288,
			-0.3536753,9.398629,
			-0.1691974,9.96429,
			2.170757,12.05057,
			1.501028,11.43888
		;
		nonmissingness.resize( X.rows(), X.cols() ) ;
		nonmissingness.setOnes() ;

		{
			Eigen::RowVectorXd expected_mean( 2 ) ;
			expected_mean << 0.1938221, 10.1426814 ;
			Eigen::MatrixXd expected_cov( 2, 2 ) ;
			expected_cov <<
				1.427126, 1.474251,
				1.474251, 1.542053
			;
			compute_mean_and_covariance( X, nonmissingness, mean, covariance ) ;
			BOOST_CHECK_SMALL( ( mean - expected_mean ).array().abs().maxCoeff(), 0.00001 ) ;
			BOOST_CHECK_SMALL( ( covariance - expected_cov ).array().abs().maxCoeff(), 0.00001 ) ;
		}

		{
			Eigen::RowVectorXd expected_mean( 20 ) ;
			expected_mean << 4.705411,5.8464,4.700252,3.382836,6.858337,5.343843,5.566398,6.788766,2.625346,5.871529,4.307139,6.023385,4.501208,5.63985,3.386614,4.817082,4.522477,4.897546,7.110664,6.469952 ;
			Eigen::MatrixXd expected_cov( 20, 20 ) ;
			expected_cov <<
				50.31921,49.39001,49.6946,49.77714,50.57331,49.45529,49.75646,50.54351,48.87431,50.48429,48.4916,50.18825,51.29384,50.72286,49.96307,49.37893,48.91692,50.82891,49.5565,49.8476,
				49.39001,48.47796,48.77694,48.85795,49.63942,48.54204,48.83765,49.61017,47.97179,49.55204,47.59615,49.26147,50.34664,49.78621,49.04045,48.46709,48.01362,49.8903,48.64139,
				48.92711,49.6946,48.77694,49.07776,49.15927,49.94556,48.84141,49.13884,49.91613,48.26764,49.85764,47.88969,49.56528,50.65714,50.09325,49.34289,48.766,48.30973,50.19799,48.94137,
				49.22885,49.77714,48.85795,49.15927,49.24091,50.02851,48.92253,49.22046,49.99903,48.34781,49.94045,47.96922,49.6476,50.74127,50.17645,49.42484,48.84699,48.38996,50.28136,
				49.02265,49.31061,50.57331,49.63942,49.94556,50.02851,50.8287,49.70503,50.00773,50.79875,49.12112,50.73923,48.73648,50.4417,51.55287,50.97901,50.21538,49.62829,49.16395,
				51.0856,49.80676,50.09933,49.45529,48.54204,48.84141,48.92253,49.70503,48.6062,48.9022,49.67574,48.03519,49.61753,47.65906,49.32658,50.41318,49.85201,49.10526,48.53115,
				48.07708,49.95624,48.70568,48.99178,49.75646,48.83765,49.13884,49.22046,50.00773,48.9022,49.20001,49.97826,48.32772,49.9197,47.9493,49.62697,50.72019,50.1556,49.40431,
				48.8267,48.36986,50.26047,49.00229,49.29013,50.54351,49.61017,49.91613,49.99903,50.79875,49.67574,49.97826,50.76882,49.09217,50.70933,48.70776,50.41197,51.52249,50.94896,
				50.18579,49.59904,49.13498,51.05549,49.77741,50.0698,48.87431,47.97179,48.26764,48.34781,49.12112,48.03519,48.32772,49.09217,47.4709,49.03465,47.09918,48.74711,49.82095,
				49.26637,48.5284,47.96103,47.51229,49.36938,48.13351,48.41624,50.48429,49.55204,49.85764,49.94045,50.73923,49.61753,49.9197,50.70933,49.03465,50.64991,48.65069,50.35291,
				51.46212,50.88927,50.12699,49.54092,49.07741,50.99567,49.71909,50.01114,48.4916,47.59615,47.88969,47.96922,48.73648,47.65906,47.9493,48.70776,47.09918,48.65069,46.73038,
				48.3654,49.43083,48.88059,48.1484,47.58547,47.14025,48.9828,47.7566,48.03712,50.18825,49.26147,49.56528,49.6476,50.4417,49.32658,49.62697,50.41197,48.74711,50.35291,
				48.3654,50.05764,51.16035,50.59086,49.83304,49.25042,48.78962,50.69664,49.42754,49.71787,51.29384,50.34664,50.65714,50.74127,51.55287,50.41318,50.72019,51.52249,49.82095,
				51.46212,49.43083,51.16035,52.28734,51.70531,50.9308,50.33534,49.86439,51.81342,50.51636,50.8131,50.72286,49.78621,50.09325,50.17645,50.97901,49.85201,50.1556,50.94896,
				49.26637,50.88927,48.88059,50.59086,51.70531,51.12975,50.36387,49.77504,49.30933,51.23666,49.95404,50.24747,49.96307,49.04045,49.34289,49.42484,50.21538,49.10526,49.40431
				,50.18579,48.5284,50.12699,48.1484,49.83304,50.9308,50.36387,49.60946,49.02945,48.57071,50.46917,49.20577,49.4948,49.37893,48.46709,48.766,48.84699,49.62829,48.53115,48.8267,
				49.59904,47.96103,49.54092,47.58547,49.25042,50.33534,49.77504,49.02945,48.45622,48.00285,49.87911,48.63048,48.91613,48.91692,48.01362,48.30973,48.38996,49.16395,48.07708,
				48.36986,49.13498,47.51229,49.07741,47.14025,48.78962,49.86439,49.30933,48.57071,48.00285,47.55372,49.41243,48.17548,48.45846,50.82891,49.8903,50.19799,50.28136,51.0856,
				49.95624,50.26047,51.05549,49.36938,50.99567,48.9828,50.69664,51.81342,51.23666,50.46917,49.87911,49.41243,51.34379,50.05849,50.35253,49.5565,48.64139,48.94137,49.02265,
				49.80676,48.70568,49.00229,49.77741,48.13351,49.71909,47.7566,49.42754,50.51636,49.95404,49.20577,48.63048,48.17548,50.05849,48.80536,49.09205,49.8476,48.92711,49.22885,
				49.31061,50.09933,48.99178,49.29013,50.0698,48.41624,50.01114,48.03712,49.71787,50.8131,50.24747,49.4948,48.91613,48.45846,50.35253,49.09205,49.38041 ;
			compute_mean_and_covariance( X.transpose(), nonmissingness.transpose(), mean, covariance ) ;
			BOOST_CHECK_SMALL( ( mean - expected_mean ).array().abs().maxCoeff(), 0.00001 ) ;
			BOOST_CHECK_SMALL( ( covariance - expected_cov ).array().abs().maxCoeff(), 0.0001 ) ;
		}
	}
}

AUTO_TEST_CASE( test_mean_and_covariance_with_no_missingness ) {
	using metro::compute_mean_and_covariance ;

	double const tolerance = 0.000000000001 ;
	Eigen::RowVectorXd mean ;
	Eigen::MatrixXd covariance ;
	Eigen::MatrixXd X ;

	X.resize( 1, 2 ) ;
	X << 1, 2 ;
	compute_mean_and_covariance( X, mean, covariance ) ;
	BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
	BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
	BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
	BOOST_CHECK_EQUAL( mean.size(), 2 ) ;
	BOOST_CHECK_EQUAL( mean(0), 1 ) ;
	BOOST_CHECK_EQUAL( mean(1), 2 ) ;
	BOOST_CHECK( std::isnan( covariance(0,0) ) ) ;
	BOOST_CHECK( std::isnan( covariance(0,1) ) ) ;
	BOOST_CHECK( std::isnan( covariance(1,0) ) ) ;
	BOOST_CHECK( std::isnan( covariance(1,1) ) ) ;

	X.resize( 2, 2 ) ;
	X <<	0, 1,
			1, 2 ;
	compute_mean_and_covariance( X, mean, covariance ) ;
	BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
	BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
	BOOST_CHECK_EQUAL( mean(0), 0.5 ) ;
	BOOST_CHECK_EQUAL( mean(1), 1.5 ) ;
	BOOST_CHECK_EQUAL( covariance(0,0), 0.5 ) ;
	BOOST_CHECK_EQUAL( covariance(0,1), 0.5 ) ;
	BOOST_CHECK_EQUAL( covariance(1,1), 0.5 ) ;
	BOOST_CHECK_EQUAL( covariance(0,1), covariance(1,0) ) ;

	X.resize( 2, 2 ) ;
	X <<	0,		1,
			0.5,	2 ;
	compute_mean_and_covariance( X, mean, covariance ) ;
	BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
	BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
	BOOST_CHECK_EQUAL( mean(0), 0.25 ) ;
	BOOST_CHECK_EQUAL( mean(1), 1.5 ) ;
	BOOST_CHECK_EQUAL( covariance(0,0), 0.125 ) ;
	BOOST_CHECK_EQUAL( covariance(0,1), 0.25 ) ;
	BOOST_CHECK_EQUAL( covariance(1,1), 0.5 ) ;
	BOOST_CHECK_EQUAL( covariance(0,1), covariance(1,0) ) ;

	compute_mean_and_covariance( X.transpose(), mean, covariance ) ;
	BOOST_CHECK_EQUAL( covariance.rows(), 2 ) ;
	BOOST_CHECK_EQUAL( covariance.cols(), 2 ) ;
	BOOST_CHECK_EQUAL( mean(0), 0.5 ) ;
	BOOST_CHECK_EQUAL( mean(1), 1.25 ) ;
	BOOST_CHECK_EQUAL( covariance(0,0), 0.5 ) ;
	BOOST_CHECK_EQUAL( covariance(0,1), 0.75 ) ;
	BOOST_CHECK_EQUAL( covariance(1,1), 1.125 ) ;
	BOOST_CHECK_EQUAL( covariance(0,1), covariance(1,0) ) ;

	X.resize( 20, 2 ) ;
	X <<
		-0.3105243,9.721345,0.9230904,10.76971,-0.2534216,9.653925,-1.579065,8.344736,1.817072,11.8996,0.4140253,10.27366,0.606559,10.52624,1.750472,11.82706,-2.246558,7.49725,0.8391379,10.90392,-0.5266163,9.140894,1.020504,11.02627,-0.6118804,9.614296,0.5836782,10.69602,-1.59382,8.367049,-0.1051234,9.739288,-0.3536753,9.398629,-0.1691974,9.96429,2.170757,12.05057,1.501028,11.43888 ;
	;

	{
		Eigen::RowVectorXd expected_mean( 2 ) ;
		expected_mean << 0.1938221, 10.1426814 ;
		Eigen::MatrixXd expected_cov( 2, 2 ) ;
		expected_cov <<
			1.427126, 1.474251,
			1.474251, 1.542053
		;
		compute_mean_and_covariance( X, mean, covariance ) ;
		BOOST_CHECK_SMALL( ( mean - expected_mean ).array().abs().maxCoeff(), 0.00001 ) ;
		BOOST_CHECK_SMALL( ( covariance - expected_cov ).array().abs().maxCoeff(), 0.00001 ) ;
	}

	{
		Eigen::RowVectorXd expected_mean( 20 ) ;
		expected_mean << 4.705411,5.8464,4.700252,3.382836,6.858337,5.343843,5.566398,6.788766,2.625346,5.871529,4.307139,6.023385,4.501208,5.63985,3.386614,4.817082,4.522477,4.897546,7.110664,6.469952 ;
		Eigen::MatrixXd expected_cov( 20, 20 ) ;
		expected_cov <<
			50.31921,49.39001,49.6946,49.77714,50.57331,49.45529,49.75646,50.54351,48.87431,50.48429,48.4916,50.18825,51.29384,50.72286,49.96307,49.37893,48.91692,50.82891,49.5565,49.8476,
			49.39001,48.47796,48.77694,48.85795,49.63942,48.54204,48.83765,49.61017,47.97179,49.55204,47.59615,49.26147,50.34664,49.78621,49.04045,48.46709,48.01362,49.8903,48.64139,
			48.92711,49.6946,48.77694,49.07776,49.15927,49.94556,48.84141,49.13884,49.91613,48.26764,49.85764,47.88969,49.56528,50.65714,50.09325,49.34289,48.766,48.30973,50.19799,48.94137,
			49.22885,49.77714,48.85795,49.15927,49.24091,50.02851,48.92253,49.22046,49.99903,48.34781,49.94045,47.96922,49.6476,50.74127,50.17645,49.42484,48.84699,48.38996,50.28136,
			49.02265,49.31061,50.57331,49.63942,49.94556,50.02851,50.8287,49.70503,50.00773,50.79875,49.12112,50.73923,48.73648,50.4417,51.55287,50.97901,50.21538,49.62829,49.16395,
			51.0856,49.80676,50.09933,49.45529,48.54204,48.84141,48.92253,49.70503,48.6062,48.9022,49.67574,48.03519,49.61753,47.65906,49.32658,50.41318,49.85201,49.10526,48.53115,
			48.07708,49.95624,48.70568,48.99178,49.75646,48.83765,49.13884,49.22046,50.00773,48.9022,49.20001,49.97826,48.32772,49.9197,47.9493,49.62697,50.72019,50.1556,49.40431,
			48.8267,48.36986,50.26047,49.00229,49.29013,50.54351,49.61017,49.91613,49.99903,50.79875,49.67574,49.97826,50.76882,49.09217,50.70933,48.70776,50.41197,51.52249,50.94896,
			50.18579,49.59904,49.13498,51.05549,49.77741,50.0698,48.87431,47.97179,48.26764,48.34781,49.12112,48.03519,48.32772,49.09217,47.4709,49.03465,47.09918,48.74711,49.82095,
			49.26637,48.5284,47.96103,47.51229,49.36938,48.13351,48.41624,50.48429,49.55204,49.85764,49.94045,50.73923,49.61753,49.9197,50.70933,49.03465,50.64991,48.65069,50.35291,
			51.46212,50.88927,50.12699,49.54092,49.07741,50.99567,49.71909,50.01114,48.4916,47.59615,47.88969,47.96922,48.73648,47.65906,47.9493,48.70776,47.09918,48.65069,46.73038,
			48.3654,49.43083,48.88059,48.1484,47.58547,47.14025,48.9828,47.7566,48.03712,50.18825,49.26147,49.56528,49.6476,50.4417,49.32658,49.62697,50.41197,48.74711,50.35291,
			48.3654,50.05764,51.16035,50.59086,49.83304,49.25042,48.78962,50.69664,49.42754,49.71787,51.29384,50.34664,50.65714,50.74127,51.55287,50.41318,50.72019,51.52249,49.82095,
			51.46212,49.43083,51.16035,52.28734,51.70531,50.9308,50.33534,49.86439,51.81342,50.51636,50.8131,50.72286,49.78621,50.09325,50.17645,50.97901,49.85201,50.1556,50.94896,
			49.26637,50.88927,48.88059,50.59086,51.70531,51.12975,50.36387,49.77504,49.30933,51.23666,49.95404,50.24747,49.96307,49.04045,49.34289,49.42484,50.21538,49.10526,49.40431
			,50.18579,48.5284,50.12699,48.1484,49.83304,50.9308,50.36387,49.60946,49.02945,48.57071,50.46917,49.20577,49.4948,49.37893,48.46709,48.766,48.84699,49.62829,48.53115,48.8267,
			49.59904,47.96103,49.54092,47.58547,49.25042,50.33534,49.77504,49.02945,48.45622,48.00285,49.87911,48.63048,48.91613,48.91692,48.01362,48.30973,48.38996,49.16395,48.07708,
			48.36986,49.13498,47.51229,49.07741,47.14025,48.78962,49.86439,49.30933,48.57071,48.00285,47.55372,49.41243,48.17548,48.45846,50.82891,49.8903,50.19799,50.28136,51.0856,
			49.95624,50.26047,51.05549,49.36938,50.99567,48.9828,50.69664,51.81342,51.23666,50.46917,49.87911,49.41243,51.34379,50.05849,50.35253,49.5565,48.64139,48.94137,49.02265,
			49.80676,48.70568,49.00229,49.77741,48.13351,49.71909,47.7566,49.42754,50.51636,49.95404,49.20577,48.63048,48.17548,50.05849,48.80536,49.09205,49.8476,48.92711,49.22885,
			49.31061,50.09933,48.99178,49.29013,50.0698,48.41624,50.01114,48.03712,49.71787,50.8131,50.24747,49.4948,48.91613,48.45846,50.35253,49.09205,49.38041 ;
		compute_mean_and_covariance( X.transpose(), mean, covariance ) ;
		BOOST_CHECK_SMALL( ( mean - expected_mean ).array().abs().maxCoeff(), 0.00001 ) ;
		BOOST_CHECK_SMALL( ( covariance - expected_cov ).array().abs().maxCoeff(), 0.0001 ) ;
	}
}
