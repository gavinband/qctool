
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include "genfile/SNPIdentifyingData.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "snptest/PerSnpFrequentistTest.hpp"
#include "snptest/NullPerSnpFrequentistTest.hpp"
#include "snptest/case_control/FrequentistTest.hpp"

namespace snptest {
	PerSnpFrequentistTest::UniquePtr PerSnpFrequentistTest::create(
		genfile::SNPIdentifyingData const& snp,
		appcontext::OptionProcessor const& options
	) {
		PerSnpFrequentistTest::UniquePtr result ;
		if( !snp.get_position().chromosome().is_sex_determining() ) {
			result.reset(
				new case_control::FrequentistTest(
					0.1,
					options.check_if_option_was_supplied( "-mimic-snptest" )
				)
			) ;
		}
		else {
			result.reset( new NullPerSnpFrequentistTest() ) ;
		}
		return result ;
	}
	
	PerSnpFrequentistTest::Results::Results():
		null_loglikelihood( std::numeric_limits< double >::quiet_NaN() ),
		alternative_loglikelihood( std::numeric_limits< double >::quiet_NaN() ),
		test_statistic( std::numeric_limits< double >::quiet_NaN() ),
		p_value( std::numeric_limits< double >::quiet_NaN() ),
		beta( std::numeric_limits< double >::quiet_NaN() ),
		standard_error( std::numeric_limits< double >::quiet_NaN() ),
		variance_covariance( Matrix::Constant( 2, 2, std::numeric_limits< double >::quiet_NaN() ) )
	{}

	PerSnpFrequentistTest::Results::Results( Results const& other ):
		null_loglikelihood( other.null_loglikelihood ),
		alternative_loglikelihood( alternative_loglikelihood ),
		test_statistic( other.test_statistic ),
		p_value( other.p_value ),
		beta( other.beta ),
		standard_error( other.standard_error ),
		variance_covariance( other.variance_covariance )
	{}

	PerSnpFrequentistTest::Results& PerSnpFrequentistTest::Results::operator=( Results const& other ) {
		null_loglikelihood = other.null_loglikelihood ;
		alternative_loglikelihood = alternative_loglikelihood ;
		test_statistic = other.test_statistic ;
		p_value = other.p_value ;
		beta = other.beta ;
		standard_error = other.standard_error ;
		variance_covariance = other.variance_covariance ;
		return *this ;
	}
}
