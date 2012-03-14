#ifndef QCTOOL_SAMPLE_SUMMARY_COMPONENT_HPP
#define QCTOOL_SAMPLE_SUMMARY_COMPONENT_HPP

#include <memory>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"

struct SampleSummaryComponent {
	static void declare_options( appcontext::OptionProcessor& options ) ;
	typedef std::auto_ptr< SampleSummaryComponent > UNiquePtr ;
	static UniquePtr create( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) ;
	SampleSummaryComponent( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) ;
	void setup( genfile::SNPDataSourceProcessor ) ;
private:
	appcontext::OptionProcessor& m_options ;
	genfile::CohortIndividualSource const& samples ;
	appcontext::UIContext& m_ui_context ;
} ;

struct SampleSummaryComputation: public boost::noncopyable {
	virtual ~SampleSummaryComputation() {}
	virtual void accumulate( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) = 0 ;
	typedef boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;
	virtual void compute( ResultCallback ) ;
} ;

struct SampleSummaryComputationManager: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
	typedef std::auto_ptr< SampleSummaryComputationManager > UniquePtr ;
	void add( std::string const& name, SampleSummaryComputation const& ) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

	void add_computation( std::string const& name, SampleSummaryComputation::UniquePtr computation ) ;

	typedef boost::function< void ( std::size_t index, genfile::SNPIdentifyingData const& snp, std::string const& computation_name, std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;
	void add_result_callback( ResultCallback ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) ;
	private:
		typedef boost::ptr_map< std::string, SampleSummaryComputation > Computations ;
		Computations m_computations ;
		std::size_t m_snp_index ;
		SNPSummaryComputation::Genotypes m_genotypes ;
		typedef boost::signals2::signal< void ( std::size_t index, genfile::SNPIdentifyingData const& snp, std::string const& computation_name, std::string const& value_name, genfile::VariantEntry const& value ) > ResultSignal ;
		ResultSignal m_result_signal ;
} ;

#endif
