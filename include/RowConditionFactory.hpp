
#ifndef __GTOOL_ROWCONDITIONFACTORY_HPP__
#define __GTOOL_ROWCONDITIONFACTORY_HPP__

#include <utility>
#include <string>
#include "RowCondition.hpp"

// exception thrown by row_condition_factory.
struct ConditionSpecException: public ConditionException
{
	ConditionSpecException( std::string const& msg ) ;
} ;

struct RowConditionFactory
{
	static std::auto_ptr< RowCondition > create_condition( std::string condition_spec ) ;

	private:
		static std::auto_ptr< RowCondition > try_to_create_compound_condition( std::string condition_spec ) ;
		static std::auto_ptr< RowCondition > try_to_create_inverted_condition( std::string condition_spec ) ;
		static std::auto_ptr< RowCondition > try_to_create_range_condition( std::string condition_spec ) ;
		static std::auto_ptr< RowCondition > create_inclusive_range_condition( std::string const& name, std::string range_spec ) ;
		static std::auto_ptr< RowCondition > create_exclusive_range_condition( std::string const& name, std::string range_spec ) ;
		static std::pair< double, double > parse_bounds( std::string const& range_spec ) ;
		} ;


#endif

