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

			double null_loglikelihood ;
			double alternative_loglikelihood ;
			double test_statistic ;
			double p_value ;
			double beta ;
			double standard_error ;
			Matrix variance_covariance ;
		} ;

	public:	
		virtual ~PerSnpFrequentistTest() {}
		virtual Results test(
			genfile::SNPIdentifyingData const& snp,
			Vector const& phenotype_values,
			FinitelySupportedFunctionSet const& all_genotypes,
			Matrix const& covariate_values,
			std::vector< std::size_t > const& indices_of_samples_to_exclude = std::vector< std::size_t >()
		) const = 0 ;
	} ;
	
}

#endif
