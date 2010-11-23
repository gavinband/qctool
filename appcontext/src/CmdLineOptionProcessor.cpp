#include <exception>
#include <string>

#include "appcontext/ProgramFlow.hpp"
#include "appcontext/CmdLineOptionProcessor.hpp"

namespace appcontext {
	CmdLineOptionProcessor::~CmdLineOptionProcessor() {
	}

	void CmdLineOptionProcessor::process( int argc, char** argv ) {
		declare_options( *this ) ;
		OptionProcessor::process( argc, argv ) ;
	}

	void CmdLineOptionProcessor::declare_options( OptionProcessor& options ) {
	}
}
