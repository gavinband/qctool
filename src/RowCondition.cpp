#include <string>
#include <sstream>
#include <cassert>

#include "GenRow.hpp"
#include "Condition.hpp"
#include "RowCondition.hpp"

void TrivialRowCondition:: format_to_stream( std::ostream& oStream ) const {
	oStream << "(1 = 1)" ;
}
