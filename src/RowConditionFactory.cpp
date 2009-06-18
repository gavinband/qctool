#include <cassert>
#include <utility>
#include <string>
#include <vector>
#include <sstream>
#include "RowConditionFactory.hpp"
#include "RowCondition.hpp"
#include "GToolException.hpp"
#include "string_utils.hpp"

ConditionSpecException::ConditionSpecException( std::string const& msg )
: ConditionException( msg )
{}


std::auto_ptr< RowCondition > RowConditionFactory::create_condition( std::string condition_spec ) {
	strip( condition_spec ) ;
	if( condition_spec.size() == 0 ) {
		return std::auto_ptr< RowCondition >( new TrivialRowCondition ) ;
	}
	
	// std::cout << "condition_factory(): " << condition_spec << ".\n" ;
	std::auto_ptr< RowCondition > condition ;

	condition = RowConditionFactory::try_to_create_compound_condition( condition_spec ) ;
	if( condition.get() != 0 )
		return condition ;
	condition = RowConditionFactory::try_to_create_inverted_condition( condition_spec ) ;
	if( condition.get() != 0 )
		return condition ;

	// If we get here, condition was not modified by AND or OR or NOT
	condition = RowConditionFactory::try_to_create_range_condition( condition_spec ) ;

    if( condition.get() == 0 ) {
        throw ConditionSpecException( "Unable to construct condition with spec \"" + condition_spec + "\"." ) ;
    }
	
	return condition ;
}



// Split up the condition by AND and OR (And having priority).
std::auto_ptr< RowCondition > RowConditionFactory::try_to_create_compound_condition( std::string condition_spec ) {
	std::vector< std::string > subconditions = split( condition_spec, "&&" ) ;
	std::auto_ptr< CompoundRowCondition > compound_condition ;

	if( subconditions.size() > 1 ) {
		compound_condition = std::auto_ptr< CompoundRowCondition >( new AndRowCondition() ) ;
		for( std::size_t i = 0 ; i < subconditions.size(); ++i ) {
			compound_condition->add_subcondition( RowConditionFactory::create_condition( strip( subconditions[i] ))) ;
		}
	}
	else {
		subconditions = split( condition_spec, "||" ) ;
		if( subconditions.size() > 1 ) {
			compound_condition = std::auto_ptr< CompoundRowCondition >( new OrRowCondition() ) ;
			for( std::size_t i = 0 ; i < subconditions.size(); ++i ) {
				compound_condition->add_subcondition( RowConditionFactory::create_condition( strip( subconditions[i] ))) ;
			}
		}
	}

	std::auto_ptr< RowCondition > condition( compound_condition.release() );
	return condition ;
}


std::auto_ptr< RowCondition > RowConditionFactory::try_to_create_inverted_condition( std::string condition_spec ) {
	std::auto_ptr< RowCondition > condition ;
	
	if( condition_spec[0] == '!' ) {
		condition_spec = condition_spec.substr( 1, condition_spec.size() - 1 ) ;
		condition = RowConditionFactory::create_condition( condition_spec ) ;
		condition = std::auto_ptr<RowCondition>( new NotRowCondition( condition )) ;
	}
	
	return condition ;
}

std::auto_ptr< RowCondition > RowConditionFactory::try_to_create_range_condition( std::string condition_spec ) {
	std::auto_ptr< RowCondition > condition ;

	std::vector<std::string> substrings = split( condition_spec, " in " ) ;
	if( substrings.size() < 2 ) {
		// range not found
		return condition ;
	}
	if( substrings.size() > 2 ) {
		throw ConditionSpecException( "Unable to construct condition with spec \"" + condition_spec + "\": spec is malformed." ) ;
	}

	// parse range.
	std::string range_str = strip( substrings[1] );
	if( range_str.size() < 2 ) {
		throw ConditionSpecException( "Unable to construct condition with spec \"" + condition_spec + "\": spec is malformed." ) ;
	}

	char first_char = range_str[0], last_char = *(range_str.rbegin()) ;
	switch( first_char ) {
		case '(':
			if( last_char != ')' ) {
				throw ConditionSpecException( "Unable to construct condition with spec \"" + condition_spec + "\": spec is malformed." ) ;
			}
			range_str = range_str.substr( 1, range_str.size() - 2 ) ;
			condition = RowConditionFactory::create_exclusive_range_condition( strip( substrings[0] ), range_str ) ;
			break ;
		case '[':
			if( last_char != ']' ) {
				throw ConditionSpecException( "Unable to construct condition with spec \"" + condition_spec + "\": spec is malformed." ) ;
			}
			range_str = range_str.substr( 1, range_str.size() - 2 ) ;
			condition = RowConditionFactory::create_inclusive_range_condition( strip( substrings[0] ), range_str ) ;
			break ;
		default:
			throw ConditionSpecException( "Unable to construct condition with spec \"" + condition_spec + "\": spec is malformed." ) ;
			break ;
	}

	return condition ;
}


std::auto_ptr< RowCondition > RowConditionFactory::create_inclusive_range_condition( std::string const& name, std::string range_spec ) {
	std::pair< double, double > bounds = RowConditionFactory::parse_bounds( range_spec ) ;
	return std::auto_ptr< RowCondition > ( new GenotypeAssayStatisticInInclusiveRange( name, bounds.first, bounds.second )) ;
}

std::auto_ptr< RowCondition > RowConditionFactory::create_exclusive_range_condition( std::string const& name, std::string range_spec ) {
	std::pair< double, double > bounds = RowConditionFactory::parse_bounds( range_spec ) ;
	return std::auto_ptr< RowCondition > ( new GenotypeAssayStatisticInExclusiveRange( name, bounds.first, bounds.second )) ;
}

// Parse bounds from the given string, specified as a pair of numbers seperated by a comma.
std::pair< double, double > RowConditionFactory::parse_bounds( std::string const& range_spec ) {
	std::vector< std::string > bounds = split( range_spec, "," ) ;
	if( bounds.size() != 2 ) {
		throw ConditionSpecException( "Unable to parse bounds with spec \"" + range_spec + "\": spec is malformed." ) ;
	}
	bounds[0] = strip( bounds[0] ) ;
	bounds[1] = strip( bounds[1] ) ;
	
	// parse the bounds
	double lower_bound, upper_bound ;
	
	

	if( bounds[0] == "-infinity" ) {
		lower_bound = -std::numeric_limits< double >::max() ;
	}
	else {
		std::istringstream aStream( bounds[0] ) ;
		aStream >> lower_bound ;
	}

	if( bounds[1] == "infinity" ) {
		upper_bound = std::numeric_limits< double >::max() ;
	}
	else {
		std::istringstream aStream( bounds[1] ) ;
		aStream >> upper_bound ;
	}

	return std::make_pair( lower_bound, upper_bound ) ;
}

