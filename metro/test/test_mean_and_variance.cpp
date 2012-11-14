
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <Eigen/Core>
#include "test_case.hpp"
#include "metro/mean_and_variance.hpp"

AUTO_TEST_CASE( test_mean_and_variance ) {
	using metro::compute_mean_and_variance ;

	double const tolerance = 0.000000000001 ;
	std::pair< double, double > result ;
	Eigen::VectorXd v ;

	v.resize( 1 ) ;
	v << 1 ;
	result = compute_mean_and_variance( v ) ;
	BOOST_CHECK_EQUAL( result.first, 1 ) ;
	BOOST_CHECK( result.second != result.second ) ;

	v.resize( 2 ) ;
	v << 1, 2 ;
	result = compute_mean_and_variance( v ) ;
	BOOST_CHECK_EQUAL( result.first, 1.5 ) ;
	BOOST_CHECK_EQUAL( result.second, 0.5 ) ;

	v.resize( 100 ) ;
	v <<
		13.869383895774,13.8555790824775,15.5394267725615,16.7007275629156,
		15.6482957546795,15.9171190143815,18.5312290051884,17.191829748269,
		15.6676974485541,15.2613314157953,17.5989588609689,14.2541697066589,
		16.7712535421831,18.809701674205,17.856615525994,17.6231877322221,
		18.0349181157301,17.9233962921262,16.9627334040908,18.6259380539988,
		14.6931989525162,16.2625662386456,16.5697048629146,15.0171721637661,
		18.4995318458278,13.0356818215237,17.4720159599778,16.9654830334714,
		18.7520337094039,21.9059280321425,16.640433064946,16.6716978031273,
		18.8343239660585,22.0483439914861,19.2359608534752,19.1794214803293,
		16.6200928091526,16.8501830141143,16.0697766432786,16.9879711617873,
		14.9729519981797,16.1861799712576,17.0509974991028,14.3082527479713,
		18.6847588999307,16.8829382772253,16.2839149850631,17.8237228829344,
		18.1833235622627,18.1823922267266,17.3308468091007,15.8865096942677,
		15.7037084964168,18.2096754180048,18.3650135484881,15.9331906934603,
		17.8280430302088,16.793392339172,16.1577061783886,16.2707888985628,
		15.7105955787009,16.2851655164224,17.2626884270657,18.2479252092362,
		18.7015305407618,15.7300702685376,16.2605452345168,17.6915504514027,
		18.4829111260908,15.5515722040466,15.899529218908,19.3905284487477,
		20.8027037805008,16.9339375404102,15.8276544734834,17.5837528570257,
		16.4857395560363,16.3764385698557,18.3459647874313,16.6776215210494,
		19.1552900524461,16.8390762478408,18.8563961389068,16.8736666126704,
		20.2330010352537,16.1236463112621,16.517233665595,16.4382610485058,
		15.2273758660282,17.3866912184788,18.7825823896346,22.1180042449853,
		18.865501747874,15.1201939099617,18.3487708924486,15.7317731422163,
		16.3978706078664,18.3442205864219,19.742176281335,17.5919706420966
	;
	
	result = compute_mean_and_variance( v ) ;
	BOOST_CHECK_CLOSE( result.first, 17.169354501255, tolerance ) ;
	BOOST_CHECK_CLOSE( result.second, 2.84840053650842, tolerance ) ;
}

AUTO_TEST_CASE( test_OnlineElementwise_mean_and_variance_no_missingness ) {
	using metro::OnlineElementwiseMeanAndVariance ;

	double const tolerance = 0.00001 ;
	std::pair< double, double > result ;

	{
		Eigen::RowVectorXd v( 1 ) ;
		v << 1 ;
		OnlineElementwiseMeanAndVariance om ;
		om.accumulate( v ) ;
		BOOST_CHECK_EQUAL( om.get_mean()(0), 1 ) ;
		BOOST_CHECK( om.get_variance()(0) != om.get_variance()(0) ) ;
		v << 2 ;
		om.accumulate( v ) ;
		BOOST_CHECK_EQUAL( om.get_mean().rows(), 1 ) ;
		BOOST_CHECK_EQUAL( om.get_mean().cols(), 1 ) ;
		BOOST_CHECK_EQUAL( om.get_variance().rows(), 1 ) ;
		BOOST_CHECK_EQUAL( om.get_variance().cols(), 1 ) ;
		BOOST_CHECK_EQUAL( om.get_mean()(0), 1.5 ) ;
		BOOST_CHECK( om.get_variance()(0) == 0.5 ) ;
	}

	{
		Eigen::MatrixXd m( 3, 2 ) ;
		m <<
			0, 1,
			1, 2,
			2, 4
		;
		OnlineElementwiseMeanAndVariance om ;
		
		for( int i = 0; i < m.rows(); ++i ) {
			om.accumulate( m.row( i ) ) ;
			std::pair< double, double > expected_result_col1 = metro::compute_mean_and_variance( m.col(0).head( i+1 ) ) ;
			std::pair< double, double > expected_result_col2 = metro::compute_mean_and_variance( m.col(1).head( i+1 ) ) ;
			BOOST_CHECK_EQUAL( om.get_mean().rows(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_mean().cols(), 2 ) ;
			BOOST_CHECK_EQUAL( om.get_variance().rows(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_variance().cols(), 2 ) ;
			BOOST_CHECK_CLOSE( om.get_mean()(0), expected_result_col1.first, tolerance ) ;
			BOOST_CHECK_CLOSE( om.get_mean()(1), expected_result_col2.first, tolerance ) ;
			if( i == 0 ) {
				BOOST_CHECK( std::isnan( om.get_variance()(0) ) ) ;
				BOOST_CHECK( std::isnan( om.get_variance()(1) ) ) ;
			} else {
				BOOST_CHECK_CLOSE( om.get_variance()(0), expected_result_col1.second, tolerance ) ;
				BOOST_CHECK_CLOSE( om.get_variance()(1), expected_result_col2.second, tolerance ) ;
			}
		}
	}

	{
		Eigen::RowVectorXd v( 1 ) ;
		v.resize( 100 ) ;
		v <<
			13.869383895774,13.8555790824775,15.5394267725615,16.7007275629156,
			15.6482957546795,15.9171190143815,18.5312290051884,17.191829748269,
			15.6676974485541,15.2613314157953,17.5989588609689,14.2541697066589,
			16.7712535421831,18.809701674205,17.856615525994,17.6231877322221,
			18.0349181157301,17.9233962921262,16.9627334040908,18.6259380539988,
			14.6931989525162,16.2625662386456,16.5697048629146,15.0171721637661,
			18.4995318458278,13.0356818215237,17.4720159599778,16.9654830334714,
			18.7520337094039,21.9059280321425,16.640433064946,16.6716978031273,
			18.8343239660585,22.0483439914861,19.2359608534752,19.1794214803293,
			16.6200928091526,16.8501830141143,16.0697766432786,16.9879711617873,
			14.9729519981797,16.1861799712576,17.0509974991028,14.3082527479713,
			18.6847588999307,16.8829382772253,16.2839149850631,17.8237228829344,
			18.1833235622627,18.1823922267266,17.3308468091007,15.8865096942677,
			15.7037084964168,18.2096754180048,18.3650135484881,15.9331906934603,
			17.8280430302088,16.793392339172,16.1577061783886,16.2707888985628,
			15.7105955787009,16.2851655164224,17.2626884270657,18.2479252092362,
			18.7015305407618,15.7300702685376,16.2605452345168,17.6915504514027,
			18.4829111260908,15.5515722040466,15.899529218908,19.3905284487477,
			20.8027037805008,16.9339375404102,15.8276544734834,17.5837528570257,
			16.4857395560363,16.3764385698557,18.3459647874313,16.6776215210494,
			19.1552900524461,16.8390762478408,18.8563961389068,16.8736666126704,
			20.2330010352537,16.1236463112621,16.517233665595,16.4382610485058,
			15.2273758660282,17.3866912184788,18.7825823896346,22.1180042449853,
			18.865501747874,15.1201939099617,18.3487708924486,15.7317731422163,
			16.3978706078664,18.3442205864219,19.742176281335,17.5919706420966
		;
	
		OnlineElementwiseMeanAndVariance om ;
	
		for( int i = 0; i < v.size(); ++i ) {
			om.accumulate( v.segment( i, 1 ) ) ;
			std::pair< double, double > expected_result = metro::compute_mean_and_variance( v.head( i+1 ) ) ;
			BOOST_CHECK_EQUAL( om.get_mean().rows(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_mean().cols(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_variance().rows(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_variance().cols(), 1 ) ;
			BOOST_CHECK_CLOSE( om.get_mean()(0), expected_result.first, tolerance ) ;
			if( i == 0 ) {
				BOOST_CHECK( std::isnan( om.get_variance()(0) ) ) ;
			} else {
				BOOST_CHECK_CLOSE( om.get_variance()(0), expected_result.second, tolerance ) ;
			}
		}
	}
}

AUTO_TEST_CASE( test_OnlineElementwise_mean_and_variance_missingness ) {
	using metro::OnlineElementwiseMeanAndVariance ;

	double const tolerance = 0.000001 ;
	std::pair< double, double > result ;


	{
		Eigen::MatrixXd m( 3, 2 ) ;
		m <<
			0, 1,
			1, 2,
			2, 4
		;

		Eigen::MatrixXd nonmissingness( 3, 2 ) ;
		nonmissingness <<
			1, 0,
			1, 1,
			0, 1
		;

		OnlineElementwiseMeanAndVariance om ;
		
		for( int i = 0; i < m.rows(); ++i ) {
			om.accumulate( m.row( i ), nonmissingness.row( i ) ) ;
			std::pair< double, double > expected_result_col1 = metro::compute_mean_and_variance( m.col(0).head( i+1 ), nonmissingness.col(0).head( i+1 ) ) ;
			std::pair< double, double > expected_result_col2 = metro::compute_mean_and_variance( m.col(1).head( i+1 ), nonmissingness.col(1).head( i+1 ) ) ;
			BOOST_CHECK_EQUAL( om.get_mean().rows(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_mean().cols(), 2 ) ;
			BOOST_CHECK_EQUAL( om.get_variance().rows(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_variance().cols(), 2 ) ;
			BOOST_CHECK_CLOSE( om.get_mean()(0), expected_result_col1.first, tolerance ) ;
			if( i == 0 ) {
				BOOST_CHECK( std::isnan( om.get_mean()(1) )) ;
			} else {
				BOOST_CHECK_CLOSE( om.get_mean()(1), expected_result_col2.first, tolerance ) ;
			}
			if( i == 0 ) {
				BOOST_CHECK( std::isnan( om.get_variance()(0) )) ;
				BOOST_CHECK( std::isnan( om.get_variance()(1) )) ;
			} else if( i == 1 ) {
				BOOST_CHECK_CLOSE( om.get_variance()(0), expected_result_col1.second, tolerance ) ;
				BOOST_CHECK( std::isnan( om.get_variance()(1) )) ;
			} else {
				BOOST_CHECK_CLOSE( om.get_variance()(0), expected_result_col1.second, tolerance ) ;
				BOOST_CHECK_CLOSE( om.get_variance()(1), expected_result_col2.second, tolerance ) ;
			}
		}
	}

	{
		Eigen::RowVectorXd v( 100 ) ;
		v <<
			13.869383895774,13.8555790824775,15.5394267725615,16.7007275629156,
			15.6482957546795,15.9171190143815,18.5312290051884,17.191829748269,
			15.6676974485541,15.2613314157953,17.5989588609689,14.2541697066589,
			16.7712535421831,18.809701674205,17.856615525994,17.6231877322221,
			18.0349181157301,17.9233962921262,16.9627334040908,18.6259380539988,
			14.6931989525162,16.2625662386456,16.5697048629146,15.0171721637661,
			18.4995318458278,13.0356818215237,17.4720159599778,16.9654830334714,
			18.7520337094039,21.9059280321425,16.640433064946,16.6716978031273,
			18.8343239660585,22.0483439914861,19.2359608534752,19.1794214803293,
			16.6200928091526,16.8501830141143,16.0697766432786,16.9879711617873,
			14.9729519981797,16.1861799712576,17.0509974991028,14.3082527479713,
			18.6847588999307,16.8829382772253,16.2839149850631,17.8237228829344,
			18.1833235622627,18.1823922267266,17.3308468091007,15.8865096942677,
			15.7037084964168,18.2096754180048,18.3650135484881,15.9331906934603,
			17.8280430302088,16.793392339172,16.1577061783886,16.2707888985628,
			15.7105955787009,16.2851655164224,17.2626884270657,18.2479252092362,
			18.7015305407618,15.7300702685376,16.2605452345168,17.6915504514027,
			18.4829111260908,15.5515722040466,15.899529218908,19.3905284487477,
			20.8027037805008,16.9339375404102,15.8276544734834,17.5837528570257,
			16.4857395560363,16.3764385698557,18.3459647874313,16.6776215210494,
			19.1552900524461,16.8390762478408,18.8563961389068,16.8736666126704,
			20.2330010352537,16.1236463112621,16.517233665595,16.4382610485058,
			15.2273758660282,17.3866912184788,18.7825823896346,22.1180042449853,
			18.865501747874,15.1201939099617,18.3487708924486,15.7317731422163,
			16.3978706078664,18.3442205864219,19.742176281335,17.5919706420966
		;
		
		Eigen::RowVectorXd nonmissingness = v ;
		nonmissingness.setConstant( 1.0 ) ;
		nonmissingness( 5 ) = 0 ;
		nonmissingness( 10 ) = 0 ;
		nonmissingness( 15 ) = 0 ;
		nonmissingness( 20 ) = 0 ;
		nonmissingness( 25 ) = 0 ;
		nonmissingness( 74 ) = 0 ;
		nonmissingness( 75 ) = 0 ;
		nonmissingness( 76 ) = 0 ;
		nonmissingness( 99 ) = 0 ;
		
		OnlineElementwiseMeanAndVariance om ;

		for( int i = 0; i < v.size(); ++i ) {
			om.accumulate( v.segment( i, 1 ), nonmissingness.segment( i, 1 ) ) ;
			std::pair< double, double > expected_result = metro::compute_mean_and_variance( v.head( i+1 ), nonmissingness.head( i+1 ) ) ;
			BOOST_CHECK_EQUAL( om.get_mean().rows(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_mean().cols(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_variance().rows(), 1 ) ;
			BOOST_CHECK_EQUAL( om.get_variance().cols(), 1 ) ;
			BOOST_CHECK_CLOSE( om.get_mean()(0), expected_result.first, tolerance ) ;
			if( i == 0 ) {
				BOOST_CHECK( std::isnan( om.get_variance()(0) )) ;
			} else {
				BOOST_CHECK_CLOSE( om.get_variance()(0), expected_result.second, tolerance ) ;
			}
		}
	}
}
