#ifndef QCTOOL_KINSHIP_COEFFICIENT_COMPUTATION_HPP
#define QCTOOL_KINSHIP_COEFFICIENT_COMPUTATION_HPP

#include <vector>
#include "Eigen/Eigen"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "SampleBySampleComputation.hpp"

// Kinship coefficient computation, as in
// Powell, Visscher and Goddard, "Reconciling the analysis of IBD and IBS in complex trait studies",
// Nature Reviews Genetics (2010).
// This agrees (up to a factor of 2) with estimate (2.2) in
// Astle & Balding, "Population Structure and Cryptic Relatedness in Genetic Association Studies", Statistical Science (2009).
// Except that Powell et al give a different expression for within-individual estimates.
// See also Rousset, "Inbreeding and relatedness coefficients: what do they measure?", Heredity (2002).
struct KinshipCoefficientComputation: public SampleBySampleComputation {
	KinshipCoefficientComputation( appcontext::OptionProcessor const& options, appcontext::UIContext& ui_context ) ;
	
	void prepare(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
	) ;
	
	// Return the number of genotypes that are called and equal between the two samples.
	double operator()(
		std::size_t const sample1,
		std::size_t const sample2,
		std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
	) ;	
	
	std::string get_summary() const ;
	
private:
	double m_threshhold ;
	std::vector< double > m_allele_frequencies ;
	Eigen::MatrixXd m_result ;
} ;


#endif
