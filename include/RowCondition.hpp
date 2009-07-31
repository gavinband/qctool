
#ifndef __GTOOL_ROWCONDITION_HPP__
#define __GTOOL_ROWCONDITION_HPP__

#include <set>
#include "GenRow.hpp"
#include "string_to_value_map.hpp"
#include "Condition.hpp"
#include "string_to_value_map.hpp"

struct ConditionValueNotFoundException: public std::exception {
	char const* what() const throw() { return "ConditionValueNotFoundException" ; }
} ;

typedef Condition< string_to_value_map > RowCondition ;
typedef CompoundCondition< string_to_value_map > CompoundRowCondition ;
typedef AndCondition< string_to_value_map > AndRowCondition ;
typedef OrCondition< string_to_value_map > OrRowCondition ;
typedef InvertCondition< string_to_value_map > NotRowCondition ;
typedef InvertCondition< string_to_value_map > InvertRowCondition ;

struct TrivialRowCondition: public RowCondition
{
	bool check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const {
		return true ;
	}
	
	void format_to_stream( std::ostream& oStream ) const ;
} ;

struct StatisticInInclusiveRange: public RowCondition
{
	StatisticInInclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const ;
	
	void format_to_stream( std::ostream& oStream ) const ;

	private:
	
		std::string const m_statistic_name ;
		double m_lower_bound, m_upper_bound, m_epsilon ;
} ;

struct StatisticInExclusiveRange: public RowCondition
{
	StatisticInExclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const ;
	
	void format_to_stream( std::ostream& oStream ) const ;
	
	private:
	
		std::string const m_statistic_name ;
		double m_lower_bound, m_upper_bound, m_epsilon ;
} ;

struct StatisticGreaterThan: public RowCondition
{
	StatisticGreaterThan( std::string const& statistic_name, double lower_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const ;
	
	void format_to_stream( std::ostream& oStream ) const ;
	
	private:
	
		std::string const m_statistic_name ;
		double m_lower_bound, m_epsilon ;
} ;

struct StatisticLessThan: public RowCondition
{
	StatisticLessThan( std::string const& statistic_name, double upper_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const ;
	
	void format_to_stream( std::ostream& oStream ) const ;
	
	private:
	
		std::string const m_statistic_name ;
		double m_upper_bound, m_epsilon ;
} ;


// Factory function for conditions.
std::auto_ptr< RowCondition > condition_factory( std::string condition_spec ) ;


#endif

