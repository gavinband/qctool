
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPUTATION_MANAGER_HPP
#define QCTOOL_SNP_SUMMARY_COMPUTATION_MANAGER_HPP

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

namespace stats {
	struct SNPSummaryComputationManager: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
	public:
		typedef std::auto_ptr< SNPSummaryComputationManager > UniquePtr ;
		typedef boost::signals2::signal<
			void (
				genfile::VariantIdentifyingData const& snp,
				std::string const& value_name,
				genfile::VariantEntry const& value
			)
		> ResultSignal ;
		typedef ResultSignal::slot_type ResultCallback ;

		typedef boost::signals2::signal<
			void (
				std::size_t sample_i,
				std::string const& value_name,
				genfile::VariantEntry const& value
			)
		> PerSampleResultSignal ;
		typedef PerSampleResultSignal::slot_type PerSampleResultCallback ;

		typedef std::map< genfile::VariantEntry, std::vector< int > > StrataMembers ;

	public:
		SNPSummaryComputationManager( genfile::CohortIndividualSource const& samples, std::string const& sex_column_name ) ;
	
		void begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) ;
		void processed_snp( genfile::VariantIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
		void end_processing_snps() ;

		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) ;
	
		void add_computation( std::string const& name, stats::SNPSummaryComputation::UniquePtr computation ) ;
		void add_result_callback( ResultCallback ) ;
		void add_per_sample_result_callback( PerSampleResultCallback ) ;

		void stratify_by( StrataMembers const&, std::string const& ) ;
	
		void set_haploid_genotype_coding( int coding ) ;
	
	private:
		genfile::CohortIndividualSource const& m_samples ;
		typedef boost::ptr_map< std::string, stats::SNPSummaryComputation > Computations ;
		Computations m_computations ;

		ResultSignal m_result_signal ;
		PerSampleResultSignal m_per_sample_result_signal ;
	
		int m_haploid_coding_column ;

		std::size_t m_snp_index ;
		SNPSummaryComputation::Genotypes m_genotypes ;
		SNPSummaryComputation::Ploidy m_ploidy ;

	private:
		std::vector< char > get_sexes( genfile::CohortIndividualSource const& samples, std::string const& sex_column_name ) const ;
		std::map< char, std::vector< int > > get_samples_by_sex( std::vector< char > const& sex ) const ;

		void fix_sex_chromosome_genotypes(
			genfile::VariantIdentifyingData const& snp,
			SNPSummaryComputation::Genotypes* genotypes,
			boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > callback
		) ;

		int determine_male_coding_column(
			genfile::VariantIdentifyingData const& snp,
			SNPSummaryComputation::Genotypes const& genotypes,
			std::vector< int > const& males
		) const ;
	} ;
}

#endif
