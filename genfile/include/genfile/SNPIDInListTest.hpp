#ifndef SNPIDINLISTTEST_HPP
#define SNPIDINLISTTEST_HPP

#include <set>
#include <string>
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	struct SNPIDInListTest: public SNPIdentifyingDataTest
	{
		SNPIDInListTest( std::set< std::string > id_fields ) ;
		bool operator()(
			std::string SNPID,
			std::string,
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
