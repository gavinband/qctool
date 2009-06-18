#include "Condition.hpp"

ConditionException::ConditionException( std::string const& msg ) 
: GToolException( msg )
{
}

