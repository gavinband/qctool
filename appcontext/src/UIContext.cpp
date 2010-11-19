#include <memory>
#include <string>
#include "appcontext/OstreamTee.hpp"
#include "appcontext/progress_bar.hpp"
#include "appcontext/CmdLineUIContext.hpp"

namespace appcontext {
	ProgressContextProxy::ProgressContextProxy( UIContext& ui_context, ProgressContextImpl const& progress_context )
		: m_ui_context( &ui_context ),
	  	  m_progress_context( &progress_context )
	{}

	ProgressContextProxy::ProgressContextProxy( ProgressContextProxy const& other )
		: m_ui_context( other.m_ui_context ),
	  	  m_progress_context( other.m_progress_context )
	{
		other.m_progress_context = 0 ;
	}

	ProgressContextProxy& ProgressContextProxy::operator=( ProgressContextProxy const& other ) {
		m_ui_context = other.m_ui_context ;
		m_progress_context = other.m_progress_context ;
		other.m_progress_context = 0 ;
		return *this ;
	}

	ProgressContextProxy::~ProgressContextProxy() {
		if( m_progress_context ) {
			m_progress_context->finish() ;
			m_ui_context->remove_progress_context( m_progress_context->name() ) ;
		}
	}

	void ProgressContextProxy::notify_progress(
		std::size_t const count,
		std::size_t const total_count
	) const {
		assert( m_progress_context ) ;
		m_progress_context->notify_progress( count, total_count ) ;
	}

	void ProgressContextProxy::finish() const {
		assert( m_progress_context ) ;
		m_progress_context->finish() ;
		m_ui_context->remove_progress_context( m_progress_context->name() ) ;
		m_progress_context = 0 ;
	}
}
