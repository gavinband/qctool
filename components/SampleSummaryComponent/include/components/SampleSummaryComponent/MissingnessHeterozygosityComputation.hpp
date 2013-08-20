
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_MISSINGNESS_HETEROZYGOSITY_COMPUTATION_HPP
#define QCTOOL_MISSINGNESS_HETEROZYGOSITY_COMPUTATION_HPP

#include <Eigen/Core>
#include <boost/optional.hpp>
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"

namespace sample_stats {
	struct MissingnessHeterozygosityComputation: public SampleSummaryComputation
	{
		MissingnessHeterozygosityComputation() ;
		MissingnessHeterozygosityComputation( genfile::Chromosome ) ;
		void accumulate( genfile::SNPIdentifyingData const&, Genotypes const&, genfile::VariantDataReader& ) ;
		void compute( ResultCallback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
	private:
		boost::optional< genfile::Chromosome > const m_chromosome ;
		std::size_t m_snp_index ;
		double m_threshhold ;
		Eigen::VectorXd m_ones ;
		Eigen::VectorXd m_total_probabilities ;
		Eigen::VectorXd m_total_calls ;
		Eigen::VectorXd m_het_snps ;
		Eigen::VectorXd m_het_snp_calls ;
	} ;
}

#endif
