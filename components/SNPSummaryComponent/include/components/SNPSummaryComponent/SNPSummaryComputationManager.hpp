
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
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

struct SNPSummaryComputationManager: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
public:
	typedef std::auto_ptr< SNPSummaryComputationManager > UniquePtr ;
public:
	
	SNPSummaryComputationManager( genfile::CohortIndividualSource const& samples ) ;
	
	void add( std::string const& name, SNPSummaryComputation const& ) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

	void add_computation( std::string const& name, SNPSummaryComputation::UniquePtr computation ) ;

	typedef boost::signals2::signal< void ( genfile::SNPIdentifyingData const& snp, std::string const& value_name, genfile::VariantEntry const& value ) > ResultSignal ;
	typedef ResultSignal::slot_type ResultCallback ;
	void add_result_callback( ResultCallback ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) ;
	
	typedef std::map< genfile::VariantEntry, std::vector< int > > StrataMembers ;
	void stratify_by( StrataMembers const&, std::string const& ) ;
	
	private:
		genfile::CohortIndividualSource const& m_samples ;
		std::vector< char > m_sexes ;
		std::map< char, std::vector< int > > m_samples_by_sex ;
		
		typedef boost::ptr_map< std::string, SNPSummaryComputation > Computations ;
		Computations m_computations ;

		ResultSignal m_result_signal ;

		std::size_t m_snp_index ;
		SNPSummaryComputation::Genotypes m_genotypes ;

	private:
		std::vector< char > get_sexes( genfile::CohortIndividualSource const& samples ) const ;
		std::map< char, std::vector< int > > get_samples_by_sex( std::vector< char > const& sexes ) const ;
		void fix_sex_chromosome_genotypes( genfile::SNPIdentifyingData const& snp, SNPSummaryComputation::Genotypes& genotypes ) const ;
		int determine_male_coding_column(
			genfile::SNPIdentifyingData const& snp,
			SNPSummaryComputation::Genotypes const& genotypes,
			std::vector< int > const& males
		) const ;
		
		
} ;

#endif
