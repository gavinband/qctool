#include <memory>
#include <string>
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/Timer.hpp"
#include "appcontext/OstreamTee.hpp"
#include "appcontext/progress_bar.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/null_ostream.hpp"

namespace appcontext {
	ApplicationContext::ApplicationContext(
		std::string const& application_name,
		std::auto_ptr< OptionProcessor > options,
		int argc,
		char** argv,
		std::string const& log_option
	 )
		: m_application_name( application_name ),
		  m_options( options ),
		  m_ui_context( new appcontext::CmdLineUIContext )
	{
		m_options->process( argc, argv ) ;
		construct_logger( log_option ) ;
		write_start_banner() ;
	}

	ApplicationContext::~ApplicationContext() {
		write_end_banner() ;
	}

	void ApplicationContext::construct_logger( std::string const& log_option ) {
		get_ui_context().logger().add_stream( "screen", std::cout ) ;
		std::string filename ;
		if( options().check_if_option_is_defined( log_option ) && options().check_if_option_has_value( log_option )) {
			std::string const& filename = options().get_value< std::string > ( log_option ) ;
			if( filename != "" ) {
				get_ui_context().logger().add_stream( "log", open_file_for_output( filename )) ;
			}
			else {
				get_ui_context().logger().add_stream( "log", std::auto_ptr< std::ostream >( new null_ostream )) ;
			}
		}
		else {
			get_ui_context().logger().add_stream( "log", std::auto_ptr< std::ostream >( new null_ostream )) ;
		}
	}

	OptionProcessor& ApplicationContext::options() const { return *m_options ; }

	ApplicationContext::UIContext& ApplicationContext::get_ui_context() const {
		return *m_ui_context ;
	}

	void ApplicationContext::write_start_banner() {
		m_ui_context->logger() << "\nWelcome to " << m_application_name << "\n"
		 	<< "(C) 2009-2010 University of Oxford\n\n";
	}

	void ApplicationContext::write_end_banner() {
		m_ui_context->logger() << "\n"
			<< "Thank you for using " << m_application_name << ".\n" ;
	}
}
