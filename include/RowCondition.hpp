
#ifndef __GTOOL_ROWCONDITION_HPP__
#define __GTOOL_ROWCONDITION_HPP__

#include <set>
#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "Condition.hpp"

typedef Condition< GenRow, GenotypeAssayStatistics > RowCondition ;
typedef CompoundCondition< GenRow, GenotypeAssayStatistics > CompoundRowCondition ;
typedef AndCondition< GenRow, GenotypeAssayStatistics > AndRowCondition ;
typedef OrCondition< GenRow, GenotypeAssayStatistics > OrRowCondition ;
typedef InvertCondition< GenRow, GenotypeAssayStatistics > NotRowCondition ;

struct TrivialRowCondition: public RowCondition
{
	bool check_if_satisfied( GenRow const& genRow, GenotypeAssayStatistics const * row_genotype_statistics_ptr ) const {
		return true ;
	}
} ;

struct GenotypeAssayStatisticInInclusiveRange: public RowCondition
{
	GenotypeAssayStatisticInInclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( GenRow const& genRow, GenotypeAssayStatistics const * row_genotype_statistics_ptr ) const ;
	
	private:
	
		std::string const m_statistic_name ;
		double m_lower_bound, m_upper_bound, m_epsilon ;
} ;

struct GenotypeAssayStatisticInExclusiveRange: public RowCondition
{
	GenotypeAssayStatisticInExclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( GenRow const& genRow, GenotypeAssayStatistics const * row_genotype_statistics_ptr ) const ;
	
	private:
	
		std::string const m_statistic_name ;
		double m_lower_bound, m_upper_bound, m_epsilon ;
} ;


// Factory function for conditions.
std::auto_ptr< RowCondition > condition_factory( std::string condition_spec ) ;


#endif

