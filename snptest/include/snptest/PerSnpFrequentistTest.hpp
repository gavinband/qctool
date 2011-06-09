#ifndef SNPTEST_PER_SNP_FREQUENTIST_TEST_HPP
#define SNPTEST_PER_SNP_FREQUENTIST_TEST_HPP

#include <memory>
#include "boost/noncopyable.hpp"
#include "Eigen/Core"
#include "appcontext/OptionProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "snptest/FinitelySupportedFunctionSet.hpp"

namespace snptest {
	struct PerSnpFrequentistTest: public boost::noncopyable {
	public:
		typedef std::auto_ptr< PerSnpFrequentistTest > UniquePtr ;
		
		// Create a test appropriate for the given snp.
		static UniquePtr create( genfile::SNPIdentifyingData const& snp, appcontext::OptionProcessor const& options ) ;

		typedef Eigen::VectorXd Point ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;

		struct Results {
			Results() ;
			Results( Results const& other ) ;
			Results& operator=( Results const& other ) ;

			double test_statistic ;
			double p_value ;
			double beta ;
			double standard_error ;
			Matrix variance_covariance ;
		} ;

	public:	
		virtual ~PerSnpFrequentistTest() {}
		virtual Results test(
			Vector const& phenotype_values,
			Matrix const& covariate_values,
			genfile::SNPIdentifyingData const& snp,
			FinitelySupportedFunctionSet const& all_genotypes
		) const = 0 ;
	} ;
	
}

#endif
