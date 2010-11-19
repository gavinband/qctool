#ifndef UICONTEXT_UICONTEXT_HPP
#define UICONTEXT_UICONTEXT_HPP

#include <string>
#include "appcontext/OstreamTee.hpp"

namespace appcontext {
	struct ProgressContextImpl
	{
		virtual void notify_progress(
			std::size_t const count,
			std::size_t const total_count
		) const = 0 ;
		virtual void finish() const = 0 ;
		virtual std::string name() const = 0 ;
		virtual ~ProgressContextImpl() {} ;
	protected:
		friend struct ProgressContextProxy ;
	} ;

	struct UIContext ;

	struct ProgressContextProxy
	{
		ProgressContextProxy( UIContext&, ProgressContextImpl const& progress_context ) ;
		ProgressContextProxy( ProgressContextProxy const& other ) ;
		ProgressContextProxy& operator=( ProgressContextProxy const& other ) ;
		~ProgressContextProxy() ;

		void notify_progress(
			std::size_t const count,
			std::size_t const total_count
		) const ;

		void finish() const ;

	private:
		UIContext* m_ui_context ;
		mutable ProgressContextImpl const* m_progress_context ;
	} ;

	struct UIContext
	{
		typedef ProgressContextProxy ProgressContext ;
	
		virtual ~UIContext() {}
		virtual OstreamTee& logger() const = 0 ;
		virtual ProgressContext get_progress_context( std::string const& name = "", std::string const& type = "bar" ) = 0 ;
	private:
		friend struct ProgressContextProxy ;
		void remove_progress_context( std::string const& name ) {
			remove_progress_context_impl( name ) ;
		}
		virtual void remove_progress_context_impl( std::string const& name ) = 0 ;
	} ;
}
#endif
