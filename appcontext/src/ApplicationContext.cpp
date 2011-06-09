#include <memory>
#include <string>
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/Timer.hpp"
#include "appcontext/OstreamTee.hpp"
#include "appcontext/progress_bar.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/null_ostream.hpp"
#include "appcontext/ProgramFlow.hpp"

namespace appcontext {
	ApplicationContext::ApplicationContext(
		std::string const& application_name,
		std::auto_ptr< OptionProcessor > options,
		int argc,
		char** argv,
		std::string const& log_option
	 ):
		m_application_name( application_name ),
		m_application_version( "" ),
		m_options( options ),
		m_ui_context( new appcontext::CmdLineUIContext )
	{
		process_options( argc, argv, log_option ) ;
		write_start_banner() ;
	}

	ApplicationContext::ApplicationContext(
		std::string const& application_name,
		std::string const& application_version,
		std::auto_ptr< OptionProcessor > options,
		int argc,
		char** argv,
		std::string const& log_option
	 ):
		m_application_name( application_name ),
		m_application_version( application_version ),
		m_options( options ),
		m_ui_context( new appcontext::CmdLineUIContext )
	{
		process_options( argc, argv, log_option ) ;
		write_start_banner() ;
	}
	
	void ApplicationContext::process_options( int argc, char** argv, std::string const& log_option ) {
		try {
			get_ui_context().logger().add_stream( "screen", std::cout ) ;
			m_options->process( argc, argv ) ;
			construct_logger( log_option ) ;
		}
		catch( OptionProcessorMutuallyExclusiveOptionsSuppliedException const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "):\n" ;
			get_ui_context().logger() << "Options \"" << e.first_option()
			<< "\" and \"" << e.second_option()
			<< "\" cannot be supplied at the same time.\n"
			<< "Please consult the documentation, or use \""
			<< m_application_name << " " << m_options->get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( OptionProcessorImpliedOptionNotSuppliedException const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "):\n" ;
			get_ui_context().logger() << "When using option \"" << e.first_option()
			<< "\", option \"" << e.second_option()
			<< "\" must also be supplied.\n"
			<< "Please consult the documentation, or use \""
			<< m_application_name << " " << m_options->get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( OptionProcessingException const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): " << e.message() << ".\n" ;
			throw HaltProgramWithReturnCode( 0 );
		}
		catch( OptionProcessorHelpRequestedException const& ) {
			get_ui_context().logger() << "Usage: "
			<< m_application_name << " <options>\n"
			<< "\nOPTIONS:\n"
			<< *(m_options)
			<< "\n" ;
			throw HaltProgramWithReturnCode( 0 );
		}
		catch( FileNotOpenedError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): Could not open the file \"" + e.filename() + "\".\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( std::exception const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): \n";
			get_ui_context().logger() << "Please use \""
			<< m_application_name << " " << m_options->get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		
	}

	ApplicationContext::~ApplicationContext() {
		write_end_banner() ;
	}

	void ApplicationContext::construct_logger( std::string const& log_option ) {
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
		m_ui_context->logger() << "\nWelcome to " << m_application_name << "\n" ;
		if( m_application_version != "" ) {
			m_ui_context->logger() << "(revision: " << m_application_version << ")\n" ;
		}
		m_ui_context->logger() << "\n(C) 2009-2010 University of Oxford\n\n";
	}

	void ApplicationContext::write_end_banner() {
		m_ui_context->logger() << "\n"
			<< "Thank you for using " << m_application_name << ".\n" ;
	}
	
	std::string const& ApplicationContext::application_name() const {
		return m_application_name ;
	}
	
}
