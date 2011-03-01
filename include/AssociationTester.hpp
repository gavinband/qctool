#ifndef ASSOCIATION_TESTER_HPP
#define ASSOCIATION_TESTER_HPP

#include <boost/math/distributions/chi_squared.hpp>
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "snptest/SNPTEST2NullModel.hpp"

struct AssociationTester: public genfile::SNPDataSourceProcessor::Callback
{
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef genfile::SingleSNPGenotypeProbabilities SingleSNPGenotypeProbabilities ;
	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;

	AssociationTester(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		appcontext::UIContext& ui_context
	) ;
	
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( SNPIdentifyingData const& id_data, SingleSNPGenotypeProbabilities const& genotypes ) ;
	void end_processing_snps() ;
private:
	appcontext::OptionProcessor const& m_options ;
	genfile::CohortIndividualSource const& m_samples ;
	appcontext::UIContext& m_ui_context ;
	std::vector< std::string > m_phenotypes ;
	std::vector< std::vector< std::size_t > > m_indices_of_samples_to_include ;
	typedef snptest2::NullModelLogLikelihood::Vector Point ;
	typedef snptest2::NullModelLogLikelihood::Vector Vector ;
	typedef snptest2::NullModelLogLikelihood::Matrix Matrix ;
	std::vector< Vector > m_phenotype_values ;
	bool m_fill_null_genotypes ;

	boost::math::chi_squared_distribution< double > m_chi_squared ;

	statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
private:
	std::vector< std::string > get_phenotypes(
		genfile::CohortIndividualSource const& samples,
		std::string const& phenotype_spec
	) const ;
		
	std::vector< std::vector< std::size_t > > get_indices_of_samples_to_include(
		std::vector< std::string > const& phenotypes,
		genfile::CohortIndividualSource const& samples
	) const ;
	
	std::vector< AssociationTester::Vector > get_phenotype_values(
		std::vector< std::string > const& phenotypes,
		genfile::CohortIndividualSource const& samples,
		std::vector< std::vector< std::size_t > > indices_of_samples_to_includes
	) const ;

	struct FrequentistTestResults {
		double test_statistic ;
		double p_value ;
		double beta ;
		double standard_error ;
		Matrix variance_covariance ;
	} ;

	FrequentistTestResults get_frequentist_test_results(
		Vector const& phenotype_values,
		SingleSNPGenotypeProbabilities const&,
		std::vector< std::size_t > const& indices_of_samples_to_include
	) const ;
	
	// forbid copying
	AssociationTester( AssociationTester const& other ) ;
} ;

#endif