#ifndef SNPIDFIELDSINLISTTEST_HPP
#define SNPIDFIELDSINLISTTEST_HPP

#include <set>
#include <string>
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	struct SNPIDFieldsInListTest: public SNPIdentifyingDataTest
	{
		SNPIDFieldsInListTest( std::set< std::string > id_fields ) ;
		bool operator()(
			std::string SNPID,
			std::string RSID,
			GenomePosition,
			char,
			char
		) const ;
	private:
		std::set< std::string > m_id_fields ;
	} ;
}

#endif
