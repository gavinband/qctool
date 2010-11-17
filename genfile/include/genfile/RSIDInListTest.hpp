#ifndef RSIDINLISTTEST_HPP
#define RSIDINLISTTEST_HPP

#include <set>
#include <string>
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	struct RSIDInListTest: public SNPIdentifyingDataTest
	{
		RSIDInListTest( std::set< std::string > id_fields ) ;
		bool operator()(
			std::string,
			std::string RSID,
			GenomePosition,
			char,
			char
		) const ;
		
		std::string display() const ;
	private:
		std::set< std::string > m_id_fields ;
	} ;
}

#endif
