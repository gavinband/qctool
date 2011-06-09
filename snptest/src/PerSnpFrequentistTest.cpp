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
		test_statistic( std::numeric_limits< double >::quiet_NaN() ),
		p_value( std::numeric_limits< double >::quiet_NaN() ),
		beta( std::numeric_limits< double >::quiet_NaN() ),
		standard_error( std::numeric_limits< double >::quiet_NaN() ),
		variance_covariance( Matrix::Constant( 2, 2, std::numeric_limits< double >::quiet_NaN() ) )
	{}

	PerSnpFrequentistTest::Results::Results( Results const& other ):
		test_statistic( other.test_statistic ),
		p_value( other.p_value ),
		beta( other.beta ),
		standard_error( other.standard_error ),
		variance_covariance( other.variance_covariance )
	{}

	PerSnpFrequentistTest::Results& PerSnpFrequentistTest::Results::operator=( Results const& other ) {
		test_statistic = other.test_statistic ;
		p_value = other.p_value ;
		beta = other.beta ;
		standard_error = other.standard_error ;
		variance_covariance = other.variance_covariance ;
		return *this ;
	}
}
