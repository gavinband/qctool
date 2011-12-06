#ifndef QCTOOL_HAPLOTYPE_FREQUENCY_COMPONENT_HPP
#define QCTOOL_HAPLOTYPE_FREQUENCY_COMPONENT_HPP

#include <string>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/RandomAccessSNPDataSource.hpp"
#include "appcontext/OptionProcessor.hpp"

// Loglikelihood of table of genotypes at 2 SNPs, given parameters specifying the haplotype frequencies.
// The parameters are \pi_01, \pi_10, and \pi_11.  Then \pi_00 is one minus the sum of the others.
// \pi_ab is the frequency of the haplotype with genotype a at SNP 1 and b at SNP 2.
// Now the probability of any table is obtained
// 
struct HaplotypeFrequencyLogLikelihood {
	typedef Eigen::MatrixXd Matrix ;
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::RowVectorXd RowVector ;

	HaplotypeFrequencyLogLikelihood( Matrix const& genotype_table ) ;
	void evaluate_at( Vector const& pi ) ;
	double get_value_of_function() const ;
	Vector get_value_of_first_derivative() ;
	Matrix get_value_of_second_derivative() ;
	private:
		Matrix const& m_genotype_table ;
		std::vector< std::vector< RowVector > > m_dpi ;
		double m_ll ;
		RowVector m_D_ll ;
		Matrix m_DDt_ll ;
} ;

struct HaplotypeFrequencyComponent: public genfile::SNPDataSourceProcessor::Callback {
public:
	typedef std::auto_ptr< HaplotypeFrequencyComponent > UniquePtr ;

	static void declare_options( appcontext::OptionProcessor& options ) ;
	static UniquePtr create( appcontext::OptionProcessor const& options ) ;

public:
	HaplotypeFrequencyComponent( genfile::SNPDataSource::UniquePtr source ) ;

public:
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const& target_snp, genfile::VariantDataReader& target_data_reader ) ;
	void compute_ld_measures(
		genfile::SNPIdentifyingData const& source_snp,
		genfile::VariantDataReader& source_data_reader,
		genfile::SNPIdentifyingData const& target_snp,
		genfile::VariantDataReader& target_data_reader
	) ;
	
	void compute_ld_measures(
		genfile::SNPIdentifyingData const& source_snp,
		genfile::SNPIdentifyingData const& target_snp,
		std::vector< std::vector< int > > const& genotypes
	) ;

	void end_processing_snps() ;
	
	typedef boost::function< void( genfile::SNPIdentifyingData const& source, genfile::SNPIdentifyingData const& target, Eigen::VectorXd const& ) > ResultCallback ;
	void send_results_to( ResultCallback callback ) ;

private:
	
	genfile::SNPDataSource::UniquePtr m_source ;
	double const m_threshhold ;
	boost::signals2::signal< void( genfile::SNPIdentifyingData const& source, genfile::SNPIdentifyingData const& target, Eigen::VectorXd const& ) > m_result_signal ;
} ;

#endif
