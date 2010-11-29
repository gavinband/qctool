#ifndef GENFILE_SNPIDMATCH_TEST_HPP
#define GENFILE_SNPIDMATCH_TEST_HPP

#include <string>
#include "SNPIdentifyingDataTest.hpp"

namespace genfile {
	struct SNPIDMatchTest: public SNPIdentifyingDataTest
	{
	public:
		SNPIDMatchTest( std::string const& expression ) ;

		bool operator()(
			std::string SNPID,
			std::string RSID,
			GenomePosition position,
			char first_allele,
			char second_allele
		) const ;
		
		std::string display() const ;

	private:
		void setup( std::string const& ) ;
		bool match( std::string const& ) const ;
	private:
		enum Type { eSNPID = 1, eRSID = 2 } ;	
		Type m_type ;
		std::string m_expression, m_prefix, m_suffix ;
		char m_wildcard_char ;
	} ;
}

#endif

