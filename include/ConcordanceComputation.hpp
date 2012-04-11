
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CONCORDANCE_COMPUTATION_HPP
#define QCTOOL_CONCORDANCE_COMPUTATION_HPP

#include <vector>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "SampleBySampleComputation.hpp"

struct ConcordanceComputation: public SampleBySampleComputation {
	ConcordanceComputation( appcontext::OptionProcessor const& options, appcontext::UIContext& ui_context ) ;
	
	void prepare(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
	) {}
	
	// Return the number of genotypes that are called and equal between the two samples.
	double operator()(
		std::size_t const sample1,
		std::size_t const sample2,
		std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
	) ;	
	
	std::string get_summary() const ;
	
private:
	double m_threshhold ;
} ;

struct PairwiseNonMissingnessComputation: public SampleBySampleComputation {
	PairwiseNonMissingnessComputation( appcontext::OptionProcessor const& options, appcontext::UIContext& ui_context ) ;
	
	void prepare(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
	) {}
	
	// Return the number of genotypes that are not missing in either sample.
	double operator()(
		std::size_t const sample1,
		std::size_t const sample2,
		std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
	) ;	

	std::string get_summary() const ;

private:
	double m_threshhold ;
} ;

#endif
