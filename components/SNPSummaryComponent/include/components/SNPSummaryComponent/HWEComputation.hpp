
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_HWE_COMPUTATION_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_HWE_COMPUTATION_HPP

#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "genfile/Error.hpp"
#include "metro/likelihood/Multinomial.hpp"
#include "components/SNPSummaryComponent/SNPHWE.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

namespace snp_summary_component {
	struct HWEComputation: public SNPSummaryComputation
	{
		HWEComputation() ;
	
		void operator()( VariantIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sexes, genfile::VariantDataReader&, ResultCallback callback ) ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	private:
		void autosomal_test( VariantIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) ;
		void autosomal_exact_test( VariantIdentifyingData const& snp, Eigen::VectorXd const& genotype_counts, ResultCallback callback ) ;
		void autosomal_multinomial_test( VariantIdentifyingData const& snp, Eigen::VectorXd const& genotype_counts, ResultCallback callback ) ;
		void X_chromosome_test( VariantIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const& sexes, ResultCallback callback ) ;

	private:
		double const m_threshhold ;
		boost::math::chi_squared_distribution< double > m_chi_squared_1df ;	
		boost::math::chi_squared_distribution< double > m_chi_squared_2df ;	
	} ;	
}

#endif
