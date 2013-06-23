
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
#include "qcdb/DBOutputter.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"

namespace impl {
	
	struct SNPDataSourceIndex {
		typedef std::auto_ptr< SNPDataSourceIndex > UniquePtr ;
		virtual ~SNPDataSourceIndex() {}
		virtual void add_index_entry( genfile::SNPIdentifyingData2 const& snp, genfile::SNPDataSink& sink ) = 0 ;
		virtual void finalise() = 0 ;
	} ;

	struct SNPOutputter: public genfile::SNPDataSourceProcessor::Callback {
	public:
		typedef std::auto_ptr< SNPOutputter > UniquePtr ;
		static UniquePtr create( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink ) ;
	public:
		SNPOutputter( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink ) ;
		~SNPOutputter() ;
		
		void begin_processing_snps( std::size_t number_of_samples ) ;
		void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
		void end_processing_snps() ;
		
		typedef boost::function< void ( genfile::SNPIdentifyingData const& snp, genfile::SNPDataSink& sink ) > IndexCallback ;
		void send_index_to( impl::SNPDataSourceIndex::UniquePtr index ) ;

	private:
		genfile::CohortIndividualSource const& m_samples ;
		bool m_manage ;
		genfile::SNPDataSink* m_sink ;
		impl::SNPDataSourceIndex::UniquePtr m_index ;
		Eigen::MatrixXd m_genotypes ;
		
	} ;
}

struct SNPOutputComponent: public boost::noncopyable {
	static void declare_options( appcontext::OptionProcessor& options ) ;

	SNPOutputComponent(
		genfile::CohortIndividualSource const& samples,
		appcontext::OptionProcessor const& options,
		appcontext::UIContext& ui_context
	) ;
	
	void setup( genfile::SNPDataSink& sink, genfile::SNPDataSourceProcessor& processor ) ;

	private:
		genfile::CohortIndividualSource const& m_samples ;
		appcontext::OptionProcessor const& m_options ;
		appcontext::UIContext& m_ui_context ;
} ;

#endif
