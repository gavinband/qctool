
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_INFO_COMPUTATION_HPP
#define QCTOOL_INFO_COMPUTATION_HPP

#include <string>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantDataReader.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "metro/SampleRange.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

namespace stats {
	namespace impl {
		struct InfoComputation: public boost::noncopyable {
		public:
			typedef SNPSummaryComputation::VariantIdentifyingData VariantIdentifyingData ;
			typedef SNPSummaryComputation::Genotypes Genotypes ;
			typedef SNPSummaryComputation::Ploidy Ploidy ;

		public:
			InfoComputation() ;
			
			void compute(
				VariantIdentifyingData const& snp,
				Genotypes const& genotypes,
				Ploidy const& ploidy
			) ;

			void compute(
				VariantIdentifyingData const& snp,
				Genotypes const& genotypes,
				Ploidy const& ploidy,
				std::vector< metro::SampleRange > const& included_samples 
			) ;
		
			double info() const { return m_info ;}
			double impute_info() const { return m_impute_info ;}

		private:
			double m_info ;
			double m_impute_info ;
			Eigen::VectorXd m_diploid_fallback_distribution ;
			Eigen::VectorXd m_haploid_fallback_distribution ;
			Eigen::VectorXd const m_diploid_levels ;
			Eigen::VectorXd const m_haploid_levels ;
			
		private:
			void compute_impl(
				VariantIdentifyingData const& snp,
				Genotypes const& genotypes,
				Ploidy const& ploidy,
				std::vector< metro::SampleRange > const& included_samples 
			) ;
		} ;
	}

	struct InfoComputation: public SNPSummaryComputation {
		void operator()(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy,
			genfile::VariantDataReader&,
			ResultCallback callback
		) ;
		
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
		
	private:
		impl::InfoComputation m_computation ;
	} ;

}


#endif
