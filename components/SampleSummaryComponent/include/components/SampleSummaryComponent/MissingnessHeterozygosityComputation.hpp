
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_MISSINGNESS_HETEROZYGOSITY_COMPUTATION_HPP
#define QCTOOL_MISSINGNESS_HETEROZYGOSITY_COMPUTATION_HPP

#include <Eigen/Core>
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"

namespace sample_stats {
	struct MissingnessHeterozygosityComputation: public SampleSummaryComputation
	{
		MissingnessHeterozygosityComputation() ;
		void accumulate( genfile::SNPIdentifyingData const&, Genotypes const&, genfile::VariantDataReader& ) ;
		void compute( ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
	private:
		std::size_t m_snp_index ;
		double m_threshhold ;
		Eigen::VectorXd m_ones ;
		Eigen::VectorXd m_missing_snps ;
		Eigen::VectorXd m_missing_call_snps ;
		Eigen::VectorXd m_het_snps ;
	} ;
}

#endif
