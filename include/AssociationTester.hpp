#ifndef ASSOCIATION_TESTER_HPP
#define ASSOCIATION_TESTER_HPP

#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"

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
	statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	
	AssociationTester( AssociationTester const& other ) ;
} ;

#endif