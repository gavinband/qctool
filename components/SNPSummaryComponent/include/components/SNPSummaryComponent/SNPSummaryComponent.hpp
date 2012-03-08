#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_HPP

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "SNPSummaryComputation.hpp"

struct SNPSummaryComputationManager: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
	typedef std::auto_ptr< SNPSummaryComputationManager > UniquePtr ;
	void add( std::string const& name, SNPSummaryComputation const& ) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

	void add_computation( std::string const& name, SNPSummaryComputation::UniquePtr computation ) ;

	typedef boost::function< void ( std::size_t index, genfile::SNPIdentifyingData const& snp, std::string const& computation_name, std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;
	void add_result_callback( ResultCallback ) ;

	private:
		typedef boost::ptr_map< std::string, SNPSummaryComputation > Computations ;
		Computations m_computations ;
		std::size_t m_snp_index ;
		SNPSummaryComputation::Genotypes m_genotypes ;
		typedef boost::signals2::signal< void ( std::size_t index, genfile::SNPIdentifyingData const& snp, std::string const& computation_name, std::string const& value_name, genfile::VariantEntry const& value ) > ResultSignal ;
		ResultSignal m_result_signal ;
} ;

struct SNPSummaryComponent: public boost::noncopyable
{
public:
	static void declare_options( appcontext::OptionProcessor& options ) ;
	SNPSummaryComponent(
		genfile::CohortIndividualSource const& samples,
		appcontext::OptionProcessor const& options,
		appcontext::UIContext& ui_context
	) ;

	genfile::SNPDataSourceProcessor::Callback::UniquePtr create() const ;

private:
	SNPSummaryComputationManager::UniquePtr create_manager() const ;
	void add_computations( SNPSummaryComputationManager& manager ) const ;
	SNPSummaryComputation::UniquePtr create_computation( std::string const& name ) const ;
private:
	genfile::CohortIndividualSource const& m_samples ;
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
} ;

#endif
