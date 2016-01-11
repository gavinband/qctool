
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SAMPLE_SUMMARY_COMPUTATION_MANAGER_HPP
#define QCTOOL_SAMPLE_SUMMARY_COMPUTATION_MANAGER_HPP

#include <string>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/signals2/signal.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/Chromosome.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "components/SampleSummaryComponent/SampleStorage.hpp"

struct SampleSummaryComputationManager: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
	typedef std::auto_ptr< SampleSummaryComputationManager > UniquePtr ;
	static UniquePtr create() ;

	void add( std::string const& name, std::string const& snp_set, SampleSummaryComputation::UniquePtr ) ;

	void begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) ;
	void processed_snp( genfile::VariantIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

	void add_computation( std::string const& name, SampleSummaryComputation::UniquePtr computation ) ;

	void send_output_to( sample_stats::SampleStorage::SharedPtr outputter ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) ;

	private:
		typedef boost::ptr_map< std::pair< std::string, std::string >, SampleSummaryComputation > Computations ;
		Computations m_computations ;
		std::size_t m_snp_index ;
		SampleSummaryComputation::Genotypes m_genotypes ;
		sample_stats::SampleStorage::SharedPtr m_outputter ;

		typedef boost::signals2::signal< void (
			std::string const& computation_name,
			std::size_t sample,
			std::string const& variable,
			std::string const& description,
			genfile::VariantEntry const& value
		) > ResultSignal ;

		ResultSignal m_result_signal ;

	private:
		void add_result_callback( ResultSignal::slot_type ) ;
} ;

#endif
