
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_RELATOTRON_HPP
#define QCTOOL_RELATOTRON_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"

#include "worker/Worker.hpp"

#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "SampleBySampleComputation.hpp"

struct Relatotron: public genfile::SNPDataSourceProcessor::Callback
{
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef genfile::SingleSNPGenotypeProbabilities SingleSNPGenotypeProbabilities ;
	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;

	Relatotron( appcontext::OptionProcessor const& options, genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) ;
	
	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( SNPIdentifyingData const& id_data, genfile::VariantDataReader& genotypes ) ;
	void end_processing_snps() ;
	
	// Run the relatedness algorithm
	void process( worker::Worker* worker = 0 ) ;
	
private:
	
	appcontext::OptionProcessor const& m_options ;
	genfile::CohortIndividualSource const& m_samples ;
	appcontext::UIContext& m_ui_context ;
	std::size_t m_number_of_samples ;
	
	typedef boost::ptr_map< std::string, SampleBySampleComputation > Computations ;
	Computations m_computations ;
	typedef std::map< std::string, std::string > ComputationFiles ;
	ComputationFiles m_computation_files ;
	
	std::vector< SNPIdentifyingData > m_snps ;
	std::vector< SingleSNPGenotypeProbabilities > m_genotypes ;
	
	typedef boost::numeric::ublas::vector< double > Vector ;
	typedef boost::numeric::ublas::matrix< double > Matrix ;
	typedef boost::numeric::ublas::zero_matrix< double > ZeroMatrix ;
	typedef boost::numeric::ublas::scalar_matrix< double > ConstantMatrix ;

	std::vector< std::size_t > m_row_samples ;
	std::vector< std::size_t > m_column_samples ;
	
private:
	
	void construct_computations() ;
	
	void process_singlethreaded(
		SampleBySampleComputation& computation,
		Matrix* bf_matrix,
		std::vector< std::size_t > const& row_samples,
		std::vector< std::size_t > const& column_samples
	) ;

	void process_multithreaded(
		SampleBySampleComputation& computation,
		Matrix* bf_matrix,
		std::vector< std::size_t > const& row_samples,
		std::vector< std::size_t > const& column_samples,
		worker::Worker& worker
	) ;
	
	std::vector< std::size_t > parse_row_spec( std::string const& spec ) const ;
	
	void perform_pairwise_computations(
		SampleBySampleComputation& computation,
		Matrix* result,
		std::vector< std::size_t > const& sample1_choice,
		std::vector< std::size_t > const& sample2_choice,
		appcontext::UIContext::ProgressContext* progress_context
	) const ;
	
	void write_sample_by_sample_matrix(
		Matrix const& bf_matrix,
		std::string const& filename,
		std::string const& description,
		std::vector< std::size_t > const& row_samples,
		std::vector< std::size_t > const& column_samples
	) const ;
} ;

#endif
