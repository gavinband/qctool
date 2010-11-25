#ifndef PROGRESS_BAR_HPP
#define PROGRESS_BAR_HPP

#include <string>

namespace appcontext {
	// Given a width and a number from 0 to 1 representing the amount of
	// progress through a task, generate a progress bar representing that amount.
	std::string get_progress_bar( std::size_t width, double progress ) ;
}
#endif
