#ifndef ASSOCIATION_TESTER_HPP
#define ASSOCIATION_TESTER_HPP

#include <boost/noncopyable.hpp>
#include <boost/signals2/signal.hpp>
#include "Eigen/Core"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "snptest/FinitelySupportedFunctionSet.hpp"

struct AssociationTester: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable
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
	void processed_snp( SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

	typedef boost::function< void ( std::size_t index, genfile::SNPIdentifyingData const& snp, std::string const& computation_name, std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;
	void add_result_callback( ResultCallback ) ;

private:
	appcontext::OptionProcessor const& m_options ;
	genfile::CohortIndividualSource const& m_samples ;
	appcontext::UIContext& m_ui_context ;
	std::vector< std::string > m_phenotype_names ;

	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;

	std::vector< Vector > m_phenotypes ;
	Matrix m_covariates ;

	typedef boost::signals2::signal< void ( std::size_t index, genfile::SNPIdentifyingData const& snp, std::string const& value_name, genfile::VariantEntry const& value ) > ResultSignal ;
	ResultSignal m_result_signal ;

private:
	std::vector< std::string > get_sample_column_names(
		genfile::CohortIndividualSource const& samples,
		std::string const& spec
	) const ;
	
	std::vector< std::vector< std::size_t > > get_indices_of_samples_to_exclude(
		std::vector< std::string > const& phenotypes,
		genfile::CohortIndividualSource const& samples
	) const ;
	
	std::vector< AssociationTester::Vector > get_phenotype_values(
		std::vector< std::string > const& phenotypes,
		genfile::CohortIndividualSource const& samples
	) const ;

Vector get_genotype_levels( Matrix const& genotypes ) const ;

} ;

#endif