#ifndef SNPTEST_PER_SNP_FREQUENTIST_TEST_HPP
#define SNPTEST_PER_SNP_FREQUENTIST_TEST_HPP

#include <memory>
#include "boost/noncopyable.hpp"
#include "Eigen/Core"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"

namespace snptest {
	struct PerSnpFrequentistTest: public boost::noncopyable {
	public:
		typedef std::auto_ptr< PerSnpFrequentistTest > UniquePtr ;
		
		// Create a test appropriate for the given snp.
		static UniquePtr create( genfile::SNPIdentifyingData const& snp ) ;

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
			genfile::SingleSNPGenotypeProbabilities const& all_genotypes,
			std::vector< std::size_t > const& indices_of_samples_to_include
		) const = 0 ;
	} ;
	
}

#endif
