
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef COMPONENTS_RELATEDNESS_COMPONENT_RELATEDNESS_COMPONENT_HPP
#define COMPONENTS_RELATEDNESS_COMPONENT_RELATEDNESS_COMPONENT_HPP

#include <boost/noncopyable.hpp>
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "worker/Worker.hpp"
#include "genfile/CohortIndividualSource.hpp"

struct RelatednessComponent: public boost::noncopyable {
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef std::auto_ptr< RelatednessComponent > UniquePtr ;
	static UniquePtr create(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker,
		appcontext::UIContext& ui_context
	) ;
	
	void setup( genfile::SNPDataSourceProcessor& processor ) const ;
	
private:
	appcontext::OptionProcessor const& m_options ;
	genfile::CohortIndividualSource const& m_samples ;
	worker::Worker* m_worker ;
	appcontext::UIContext& m_ui_context ;

private:
	RelatednessComponent(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker,
		appcontext::UIContext& ui_context
	) ;
	
	std::string get_PCA_filename_prefix() const ;
} ;

#endif
