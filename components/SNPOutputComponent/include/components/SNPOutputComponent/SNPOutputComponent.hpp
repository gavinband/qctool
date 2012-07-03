
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef COMPONENTS_SNPOUTPUTCOMPONENT_SNPOUTPUTCOMPONENT_HPP
#define COMPONENTS_SNPOUTPUTCOMPONENT_SNPOUTPUTCOMPONENT_HPP

#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"

namespace impl {
	struct SNPOutputter: public genfile::SNPDataSourceProcessor::Callback {
	public:
		static UniquePtr create( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink::UniquePtr sink ) ;
		static UniquePtr create( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink ) ;
	public:
		SNPOutputter( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink::UniquePtr sink ) ;
		SNPOutputter( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink ) ;
		~SNPOutputter() ;
		
		void begin_processing_snps( std::size_t number_of_samples ) ;
		void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
		void end_processing_snps() ;
	private:
		genfile::CohortIndividualSource const& m_samples ;
		bool m_manage ;
		genfile::SNPDataSink* m_sink ;
		Eigen::MatrixXd m_genotypes ;
	} ;
}

struct SNPOutputComponent: public boost::noncopyable {
	static void setup( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink::UniquePtr sink, genfile::SNPDataSourceProcessor& processor ) ;
	static void setup( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink, genfile::SNPDataSourceProcessor& processor ) ;
} ;

#endif
