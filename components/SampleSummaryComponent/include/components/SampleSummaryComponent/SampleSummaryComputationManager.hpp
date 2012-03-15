#ifndef QCTOOL_SAMPLE_SUMMARY_COMPUTATION_MANAGER_HPP
#define QCTOOL_SAMPLE_SUMMARY_COMPUTATION_MANAGER_HPP

#include <string>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/signals2/signal.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/Chromosome.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"

struct SampleSummaryComputationManager: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
	typedef std::auto_ptr< SampleSummaryComputationManager > UniquePtr ;

	void add( std::string const& name, std::string const&, SampleSummaryComputation::UniquePtr ) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

	void add_computation( std::string const& name, SampleSummaryComputation::UniquePtr computation ) ;

	typedef boost::signals2::signal< void ( std::string const& computation_name, std::string const& snp_set, std::string const& value_name, genfile::VariantEntry const& value ) > ResultSignal ;
	void add_result_callback( ResultSignal::slot_type ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) ;
	private:
		typedef boost::ptr_map< std::pair< std::string, std::string >, SampleSummaryComputation > Computations ;
		Computations m_computations ;
		std::size_t m_snp_index ;
		SampleSummaryComputation::Genotypes m_genotypes ;
		ResultSignal m_result_signal ;
} ;

#endif
