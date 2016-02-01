
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_HAPLOTYPE_FREQUENCY_COMPONENT_HPP
#define QCTOOL_HAPLOTYPE_FREQUENCY_COMPONENT_HPP

#include <string>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "HaplotypeFrequencyLogLikelihood.hpp"

struct HaplotypeFrequencyComponent: public genfile::SNPDataSourceProcessor::Callback {
public:
	typedef std::auto_ptr< HaplotypeFrequencyComponent > UniquePtr ;

	static void declare_options( appcontext::OptionProcessor& options ) ;
	static UniquePtr create(
		genfile::CohortIndividualSource const& samples,
		std::string const& source_sample_id_column,
		genfile::CohortIndividualSource::UniquePtr ld_samples,
		std::string const& ld_sample_id_column,
		genfile::SNPDataSource::UniquePtr source,
		appcontext::OptionProcessor const& options,
		appcontext::UIContext& ui_context
	) ;

	void set_max_distance( uint64_t distance ) ;

public:
	HaplotypeFrequencyComponent( genfile::SNPDataSource::UniquePtr source, appcontext::UIContext& ui_context ) ;

public:
	void begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) ;
	void processed_snp( genfile::VariantIdentifyingData const& target_snp, genfile::VariantDataReader& target_data_reader ) ;
	void compute_ld_measures(
		genfile::VariantIdentifyingData const& source_snp,
		genfile::VariantDataReader& source_data_reader,
		genfile::VariantIdentifyingData const& target_snp,
		genfile::VariantDataReader& target_data_reader
	) ;
	
	void compute_ld_measures(
		genfile::VariantIdentifyingData const& source_snp,
		genfile::VariantIdentifyingData const& target_snp,
		std::vector< Eigen::VectorXd > const& genotypes
	) ;

	void end_processing_snps() ;
	
	typedef boost::function< void( genfile::VariantIdentifyingData const& source, genfile::VariantIdentifyingData const& target, std::string const&, genfile::VariantEntry const& ) > ResultCallback ;
	void send_results_to( ResultCallback callback ) ;

private:
	genfile::SNPDataSource::UniquePtr m_source ;
	appcontext::UIContext& m_ui_context ;
	double const m_threshhold ;
	int64_t m_max_distance ;
	boost::signals2::signal< void( genfile::VariantIdentifyingData const& source, genfile::VariantIdentifyingData const& target, std::string const&, genfile::VariantEntry const& ) > m_result_signal ;
} ;

#endif
