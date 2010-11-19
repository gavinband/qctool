#include <exception>
#include <string>

#include "appcontext/ProgramFlow.hpp"
#include "appcontext/CmdLineOptionProcessor.hpp"

namespace appcontext {
	CmdLineOptionProcessor::~CmdLineOptionProcessor() {
	}

	void CmdLineOptionProcessor::process( int argc, char** argv ) {
		try {
			declare_options( *this ) ;
			OptionProcessor::process( argc, argv ) ;
		}
		catch( OptionProcessorMutuallyExclusiveOptionsSuppliedException const& e ) {
			std::cerr << "Options \"" << e.first_option()
			<< "\" and \"" << e.second_option()
			<< "\" cannot be supplied at the same time.\n"
			<< "Please consult the documentation, or use \""
			<< get_program_name() << " " << get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( OptionProcessorImpliedOptionNotSuppliedException const& e ) {
			std::cerr << "When using option \"" << e.first_option()
			<< "\", option \"" << e.second_option()
			<< "\" must also be supplied.\n"
			<< "Please consult the documentation, or use \""
			<< get_program_name() << " " << get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( OptionProcessorHelpRequestedException const& ) {
			std::cerr << "Usage: "
			<< get_program_name() << " <options>\n"
			<< "\nOPTIONS:\n"
			<< *this
			<< "\n" ;
			throw HaltProgramWithReturnCode( 0 );
		}
		catch( std::exception const& exception ) {
			std::cerr << "!! Error: " << exception.what() << ".\n";
			std::cerr << "Please use \""
			<< get_program_name() << " " << get_help_option_name()
			<< "\" for more information.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
	}

	void CmdLineOptionProcessor::declare_options( OptionProcessor& options ) {
	}
}
