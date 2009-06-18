#ifndef __GTOOL__SAMPLECONDITION_HPP
#define __GTOOL__SAMPLECONDITION_HPP

#include "Condition.hpp"

typedef Condition< SampleRow, GenotypeProportions > SampleCondition ;

struct SampleInList: public SampleCondition
{
	SampleInList( std::string filename ) ;
	bool check_if_satisfied( SampleRow const& genRow, GenotypeProportions const * ) const ;

	protected:
		FromFileSet< std::set< std::string > > m_sample_id_list ;
} ;




#endif
